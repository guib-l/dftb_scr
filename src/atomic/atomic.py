#!/usr/bin/env python3
import __future__

# Import standard de pythons3
import os
import sys
import copy
import time

import numpy as np
import ase
from ase import Atoms

import warnings
        
import scipy
from scipy.optimize import bisect
from scipy.linalg import eig,eigh,triu,tril

from .environnement import Environnement
from src.pydftb.gamma import Gamma

class Atomic(Environnement):
    """Atomic"""
    def __init__(self,
            image:ase.Atoms,
            Hamiltonian:np.ndarray = None,
            Overlaps:np.ndarray    = None,
            Hubbard:np.ndarray     = None,
            Charges:np.ndarray     = None,
            Occupations:np.ndarray = None,
            verbose:np.int32       = 0,
            **kwargs):
        """__init__"""
        Environnement.__init__(self, image, **kwargs)
        
        self.set_net_charges(Charges, misc=False )
        if Occupations is not None:
            self.set_population(population=Occupations)

        self.ham = Hamiltonian
        self.ovr = Overlaps
        self.hub = Hubbard

        # On définit la matrice gamma ici même. 
        self.define_gamma(ShellResolvedSCC = self.ShellResolvedSCC)
        
        self.eigenVec = None
        self.eigenVal = None
        self.density  = None

        self.previous_net_charges     = None
        self.previous_partial_charges = None

        # Variables accesoires 
        self.energies = dict()

        self.set_Wa()

    def define_gamma(self, ShellResolvedSCC=True):
        """define_gamma"""

        self.gam = dict(o=None,l=None,a=None)

        if ShellResolvedSCC:
            self.gam["o"] = self._calc_gamma_heteroNuclear(self.hub, "o")
            self.gam["l"] = self._calc_gamma_heteroNuclear(self.hub, "l")
            self.gam["a"] = self._calc_gamma_heteroNuclear(self.hub, "a")

        if ShellResolvedSCC is False:
            _hub = {}
            for key,items in self.hub.items():
                _hub.update( {key : np.array([items[0]]*3) })

            self.gam["o"] = self._calc_gamma_heteroNuclear(_hub, "o")
            self.gam["l"] = self._calc_gamma_heteroNuclear(_hub, "l")
            self.gam["a"] = self._calc_gamma_heteroNuclear(_hub, "a")

        assert self.gam["o"].shape[0]==self.ham.shape[0], \
            ValueError("\n [ERROR] : La matrice Gamma employé au calcul de \
l'hamiltonien se doit d'être de a même taille que l'hamiltonien 'zeros'. Cela \
signifie qu'il est nécéssaire d'avoir 'shelltype' correspondant au charges, \
et à l'hamiltonien. \n\t[ Hamiltonien.shape == Gamma.shape ] \
\n\t[ charges.shape     == Gamma.shape[0] ] \
\n\t[ charges.shape     == shelltype ] \n")


    def __call__(self, charges=None, temperature=1e-5):
        """__call__"""
        self.reset()
        self.diagonalize(charges=charges)
        self.set_population(temperature=temperature)

        self.set_charges()
        self.set_total_energy()

        return self.partial_charges["l"]

    def reset(self):
        """reset"""
        self.density  = None
        #self.eigenVec = None
        #self.eigenVal = None


    def diagonalize(self, charges=None, eigMethod=scipy.linalg.eigh, **kwargs)->None:
        """diagonalize"""
        self.startl("Diagonalisation")
        
        shell  = "a"
        if self.ShellResolvedSCC:
            shell  = "l"
        
        if charges is None:
            _q0   = self.partial_charges[shell]
        else: _q0 = charges

        _Gamma = self.gam[shell]
        
        assert _Gamma.shape[0] == _q0.shape[0], \
            ValueError("Gamma and charge needs to be the same size")
        
        _h1    = self._calc_Hcoulomb(_Gamma, _q0)
        #_h1    = self._calc_Hcoulomb(_Gamma, self.net_charges["a"])

        _realH = self.ham + _h1 

        self._realH = _realH
        
        _realH = tril(_realH,0)
        _realS = tril(self.ovr,0)

        _e,_c = eigMethod(_realH,_realS, **kwargs) 
        
        self.set_eigenValue(_e)
        self.set_eigenVector(_c)

        #print( _c )
        #sys.exit(0)

        self.stopl("Diagonalisation")

    def _calculate_eigenValue(self, charges=None, update=False, misc=True)->None:
        """_calculate_eigenValue"""

        shell  = "a"
        if self.ShellResolvedSCC:
            shell  = "l"

        if charges is None:
            _q0   = self.partial_charges[shell]
        else: _q0 = charges.copy()

        _Gamma = self.gam[shell]
        _h1    = self._calc_Hcoulomb(_Gamma, _q0)
        _realH = self.ham + _h1 

        p = self.population.copy()

        C = np.copy(self.get_eigenVector())
        #E = np.einsum( "ia,ij,ja->a",C.T,_realH,C )
        E = np.diag( C.T @ self._realH @ C )

        #idx = np.argsort(E)
        #E,C,p = self.sort_eigenVector(E,C,p)

        if update:
            self.set_eigenValue(E)
            self.set_eigenVector(C)
            self.set_population(population=p)
        if misc:
            return 0
        else: return E

    def sort_eigenVector(self,eigenValues, eigenVectors, pop):
        """sort_eigenVector"""
        idx = np.argsort(eigenValues)
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[idx]
        pop = pop[idx]
        return eigenValues, eigenVectors, pop


#################################################
    # Public

    # *******************************************
    # Getter

    def get_density(self)->np.ndarray:
        """get_density"""
        if self.density is None:
            self.set_density()
        return self.density

    def get_eigenVector(self,)->np.ndarray:
        """get_eigenVector"""
        return self.eigenVec

    def get_eigenValue(self,)->np.ndarray:
        """get_eigenValue"""
        return self.eigenVal
    
    def get_total_energy(self):
        """get_total_energy"""
        Esum = 0.0
        for key,item in self.energies.items():
            Esum += item
        return Esum

    def get_Wa(self,)->np.ndarray:
        """get_Wa"""
        return self.Wa

    # *******************************************
    # Setter
    def set_total_energy(self,vector=None, charges=None, gamma=None,):
        """set_total_energy"""
        self.set_bsEnergy( vector=None )
        self.set_sccEnergy( charges=charges, gamma=gamma )

    def set_Wa(self,)->None:
        """set_Wa"""
        self.Wa = self._calc_Wa()

    def set_bsEnergy(self,other=None, vector=None):
        """set_bsEnergy"""
        if other is None:
            self.energies["bs_energy"] = self._calc_energies(option="bs", eigenVector=vector)
        else:
            self.energies["bs_energy"] = other

    def set_sccEnergy(self,other=None,charges=None, gamma=None):
        """set_sccEnergy"""
        if other is None:
            self.energies["scc_energy"] = self._calc_energies( 
                option="scc", charges=charges, gamma=gamma)
        else:
            self.energies["scc_energy"] = other

    def set_eigenVector(self, other)->None:
        """set_eigenVector"""
        self.eigenVec = other

    def set_eigenValue(self, other)->None:
        """set_eigenValue"""
        self.eigenVal = other

    def set_density(self,)->None:
        """set_density"""
        self.density = self._calc_density()

    def set_charges(self,)->None:
        """set_density"""
        self.previous_net_charges     = copy.deepcopy( self.net_charges )
        self.previous_partial_charges = copy.deepcopy( self.partial_charges )

        self.set_net_charges( self._calc_charges(), misc=True )

    # Overflow method
    def set_population(self, 
            temperature = 0., 
            Nelectrons  = None,  
            energies    = None, 
            population  = None, 
            update      = True,
            **kwargs )->None:
        """set_population""" 

        if population is not None:
            self.set_option("population", population)
            return

        if energies is None:
            energies   = self.get_eigenValue().copy()
        if Nelectrons == None:
            Nelectrons = np.sum(self.net_charges["o"])

        other = self._calc_fermilvl(energies, 
                temperature, Nelectrons, **kwargs)
        if update:
            self.set_option("population", other[0])
            self.set_option("fermi_level", other[1])
        else:
            return other

    # *******************************************
    # Calc'tter

    def _calc_Wa(self, shell="o")->np.ndarray:
        """ _calc_Wa """
        self.startl("Set_Wa")
        N   = len(self.ovr)
        _Wa = np.zeros((N ,*self.ham.shape))

        for i in range(N):
            bp,bm = (i+1),i
            Sa = np.zeros(self.ovr.shape, dtype=np.float64)
            Sa[:,bm:bp] = self.ovr[:,bm:bp]
            temp  = np.zeros_like(self.ovr)
            temp  = 0.5 *(Sa + Sa.T) 
            _Wa[i] = temp

        self.stopl("Set_Wa")
        return _Wa

    def _calc_energies(self, charges=None, eigenVector=None, gamma=None, option="bs"):
        """_calc_energies"""
        self.startl("Set_Energy")
        nrj = 0.
        if option == "bs":
            if eigenVector is None: _c = self.get_eigenVector()
            else: _c = eigenVector
            nrj = np.einsum( "i,ni,nm,mi->",self.population,_c,self.ham,_c )

        elif option == "scc":
            if charges is None:
                charges = self.partial_charges[self.shelltype]
            if gamma is None: gamma = self.gam[self.shelltype]
            assert gamma.shape[0]==charges.shape[0], \
                ValueError("Size of gamma and charges doesn't match !")
            nrj = np.einsum( "i,ij,j",charges,gamma,charges )*0.5

        else:
            raise NotImplementedError("Energy calculation demand isn't set.")
        self.stopl("Set_Energy")
        return nrj

    def _calc_density(self, vector=None, population=None)->None:
        """_calc_density"""
        self.startl("Set_Density")

        if population is None:
            population     = self.population.ravel()
        if vector is None:
            vector   = self.eigenVec
        
        dens = 0.5 * np.einsum("n,an,bn",population,vector,vector)
        density = (dens + dens.T)

        self.stopl("Set_Density")
        return density

    def _calc_charges(self, _density=None)->None:
        """_calc_charges"""
        self.startl("Set_Charges")

        if _density is None:
            _density = self.get_density()

        if len(_density.shape) > 2:
            _density = np.sum(_density, axis=2)

        chg = np.einsum('...ii->...i', _density @ self.ovr)

        self.stopl("Set_Charges")
        return chg 

    def _calc_Hcoulomb(self, 
            gamma:np.ndarray=None, 
            _charges:np.ndarray=None, ) -> np.ndarray:
        """_calc_Hcoulomb"""
        self.startl("Set_Hcoulomb")

        assert gamma.shape[0] == _charges.shape[0], \
                ValueError("Gamma and Charges needs to be the same size.")
        
        Ua = gamma @ _charges

        h1 = np.zeros((Ua.shape[0],Ua.shape[0]),)
        
        for i,mu in enumerate(Ua):
            for j,nu in enumerate(Ua):
                h1[i,j] = 0.5 * ( mu + nu )

        number = self.atomic_number
        tmp    = np.zeros_like(self.ovr)

        tmp = np.zeros_like(self.ovr)

        number = self.atomic_number
        if self.ShellResolvedSCC:
            number = self.azimuthal_number
        
        ci = 0
        for i,u in enumerate(number):
            cj = 0
            for j,v in enumerate(number):
                tmp[ci:ci+int(u),cj:cj+int(v)] = h1[i,j]
                cj += int(v)
            ci += int(u)
        
        h1 = tmp * self.ovr

        self.stopl("Set_Hcoulomb")
        return h1

    def _calc_fermilvl(self, 
            energies:np.ndarray    = None,
            temperature:np.float64 = 0.,
            Nelectrons:np.float64  = 0., 
            xtol:np.float64        = 1e-14)->tuple:
        """_calc_fermilvl"""
        self.startl("Set_FermiLvl")
        
        warnings.filterwarnings('ignore')
        
        
        if temperature == 0.00:
            nelec = Nelectrons
            _occ  = np.zeros(len(energies))
            for o,value in enumerate(_occ):
                if nelec>=2:
                    _occ[o] = 2
                    nelec -= 2.
                elif nelec==0.00:pass
                else:
                    _occ[o] = nelec
                    nelec = 0.00
            _mu = energies[  list(_occ).index(0.0)  ]
            
            return _occ,_mu
        

        kbT = temperature * (ase.units.kB / ase.units.Hartree)

        fermi = lambda x,mu:2*(np.exp((x-mu)/kbT)+1)**-1

        def cost(mu):
            sumF = 0
            for _e in energies: sumF += fermi(_e,mu)
            return sumF - Nelectrons

        _mu  = bisect(cost, np.min(energies), np.max(energies),rtol=xtol, xtol=xtol, maxiter=99999)

        _occ = np.array([fermi(e,_mu) for e in energies])
        _occ = np.where(np.abs(_occ)>xtol, _occ,0) 

        self.stopl("Set_FermiLvl")
        return _occ,_mu

    def _calc_gamma_heteroNuclear(self, 
            hubbard:dict  = None, 
            shelltype:str = "o") -> np.ndarray:
        """ _calc_gamma_heteroNuclear """
        
        self.startl("Set_Gamma")

        N  = self.number
        self._set_distances_nopbc(units=ase.units.Bohr**-1) 
        self.distance_matrix += np.identity(N)

        invDistance_matrix = self.distance_matrix**-1  - np.identity(N)
        Distance_matrix    = self.distance_matrix      - np.identity(N)


        # *************************************************************** #
        # Horrible boucle non optimisée 
        lst_atoms  = self.name
        it         = 1
        HubbardLst = []
        AtomsLst   = []
        if shelltype=="o":N=[1,3,5]
        if shelltype=="l":N=[1,1,1]
        if shelltype=="a":N=[1,0,0]

        for l,iat in enumerate(self.angular_momentum):

            if iat>=0:
                HubbardLst.append( hubbard[lst_atoms[l]][0] )
                AtomsLst.append(l)
            if iat>=1:
                for j in range(N[1]):
                    HubbardLst.append( hubbard[lst_atoms[l]][1] )
                    AtomsLst.append(l)
            if iat>=2:
                for j in range(N[2]):
                    HubbardLst.append( hubbard[lst_atoms[l]][2] )
                    AtomsLst.append(l)
        
        
        gamma = np.zeros((len(HubbardLst),len(HubbardLst),), dtype=np.float64) 
        for k in range(len(HubbardLst)):
            for l in range(len(HubbardLst)):
                if Distance_matrix[AtomsLst[k], AtomsLst[l]]<= 0.:
                    value = - Gamma( Distance_matrix[AtomsLst[k], AtomsLst[l]], 
                                             HubbardLst[k], HubbardLst[l])
                else:
                    value = invDistance_matrix[AtomsLst[k],AtomsLst[l]] - \
                        Gamma( Distance_matrix[AtomsLst[k], AtomsLst[l]], 
                                             HubbardLst[k], HubbardLst[l])
                #print( k,l,k%2,l%2, Gamma( Distance_matrix[AtomsLst[k], AtomsLst[l]], 
                #                             HubbardLst[k], HubbardLst[l]) )
                gamma[k,l] = value
        
        self.stopl("Set_Gamma")
        return gamma

    def _calc_hubbard( self,
            hubbard:dict  = None,
            shelltype:str = "o"):
        """_calc_hubbard"""
        HubbardLst = []
        AtomsLst   = []
        if shelltype=="o":N=[1,3,5]
        if shelltype=="l":N=[1,1,1]
        if shelltype=="a":N=[1,0,0]

        lst_atoms  = self.name
        for l,iat in enumerate(self.angular_momentum):

            if iat>=0:
                HubbardLst.append( hubbard[lst_atoms[l]][0] )
                AtomsLst.append(l)
            if iat>=1:
                for j in range(N[1]):
                    HubbardLst.append( hubbard[lst_atoms[l]][1] )
                    AtomsLst.append(l)
            if iat>=2:
                for j in range(N[1]):
                    HubbardLst.append( hubbard[lst_atoms[l]][2] )
                    AtomsLst.append(l)
        return np.array(HubbardLst,dtype=np.float64)

    def _calc_forces(self,):
        """_calc_forces"""
        return 0


    def write_atomic(self,extend_format=False, label="dftb"):
        """write_atomic"""
        
        def reshape(q,n):
            nq=[]
            c=np.c_[q,n]
            for i in range(self.number):
                nq.append(list( q[n==i] ))
            return nq
            
        properties = []
        exclude    = ["_results","timer_id","_start","_stop","_level","is_running","Wa",
                      "gam","hub","previous_partial_charges","previous_net_charges"]
        matrix     = ["distance_matrix","eigenVec","eigenVal","density",]
        if extend_format:
            matrix += ["ham","ovr",]

        try: os.remove('atomic_%s.out'%label)
        except:pass
        line = self.__dict__
        with open('atomic_%s.out'%label,"a") as fout:
            
            fout.write("*"*80+"\n")
            fout.write(" === {}\n".format("Properties of calculation"))

            special = ["net_charges", "partial_charges",]
            for key,item in line.items():
                if key in exclude:continue
                if key in matrix:continue
                if key in special:continue

                if "_valence" in key or key=="population":
                    fout.write("{:>25} : \n\t".format(key))
                    for i,a in enumerate(item):
                        fout.write("{:1.8e}  ".format(a))
                        if i%10==9:fout.write("\n\t")
                        if i==len(item)-1 and i%10!=9:fout.write("\n")
                else:
                    fout.write("{:>25} : {}\n".format(key,item))
                properties.append(key)

            for key in special:

                name = ["orbitals_name","azimuthal_name","atomic_name"]
                fout.write("{:>25} : ".format(key))
                for i,(k,it) in enumerate(line[key].items()):
                    fout.write("\n {:>5} : ".format(k))
                    tab=reshape(it,line[name[i]])
                    for g,_l in enumerate(tab):
                        for j,a in enumerate(_l):
                            fout.write("{:1.8e}  ".format(a))
                        if g!=len(tab)-1: fout.write("\n{:>9}".format(""))
                fout.write("\n")
                properties.append(key)

            
            fout.write("\n === {}\n".format("Matrix of calculation"))

            for key in ["gam",]:
                name = ["orbitals_name","azimuthal_name","atomic_name"]
                fout.write("{:>25} : ".format("Gamma"))
                for i,(k,it) in enumerate(line[key].items()):
                    fout.write("\n {:>5} : ".format(k))
                    tab = it
                    for g,_l in enumerate(tab):
                        for j,a in enumerate(_l):
                            fout.write("{:1.8e}  ".format(a))
                        if g!=len(tab)-1: fout.write("\n{:>9}".format(""))
            for key in matrix:
                fout.write("\n{:>25} : ".format(key))
                tab = line[key]
                if tab is None:
                    fout.write( "None" )
                    continue
                fout.write( "\n{:>9}".format('') )
                for g,_l in enumerate(tab):
                    if len(tab.shape) == 1:
                        fout.write("{:1.8e}  ".format(_l))
                    if len(tab.shape) == 2:
                        for j,a in enumerate(_l):
                            fout.write("{:1.8e}  ".format(a))
                        if g!=len(tab)-1: fout.write("\n{:>9}".format(""))






