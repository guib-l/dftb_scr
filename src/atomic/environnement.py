#!/usr/bin/env python3
import __future__

# Import standard de pythons3
import os
import sys
import ase
import math
import copy
import time
import numpy as np
import scipy

from src.misc.series import structure
from src.misc.timer import Timer
from src.misc.misc import (get_angular_momentum, atomic_data)

def _get_momentum(l):
    """ _get_momentum """
    if l==0:return 1
    if l==1:return 4
    if l==2:return 9

class Environnement(structure,Timer):
    """Environnement"""
    ORBITALE = ["s","py","pz","px","dxy","dyz","dz2","dxz","dx2-y2",]
    ATTRIBUT = ["name",                 # Symbole qui représente la molécule
                "number",               # Nombre atome qui compose la molécule
                "type_name",            # List des élements qui compose la 
                                        # molécule
                "type_interraction",    # Liste de toutes les interraction 
                                        # inter-espèce 
                "interraction_matrix",  # Matrice des interractions inter-espèce
                "max_angular_momentum", # 
                "net_charges",          # 
                "partial_charges",      # 
                "electron_number",      # 
                "distance_matrix",      # 
                "angular_momentum",     # 
                "atomic_number",        #  
                "atomic_name",          # 
                "azimuthal_number",     # 
                "azimuthal_name",       # 
                "orbitals_number",      # 
                "orbitals_name",        # 
                "atomic_valence",       # 
                "azimuthal_valence",    # 
                "orbitals_valence",     # 
            ]

    def __init__(self,
            image:ase.Atoms    = None,
            shelltype:str      = "o",
            ShellResolvedSCC:bool  = True,
            verbose:np.int32   = 0,
            low:bool           = False)->None:
        """__init__"""
        structure.__init__(self, image)

        Timer.__init__(self, "Calculation DFTB")

        self.image     = image
        self.verbose   = verbose
        self.shelltype = shelltype

        self.ShellResolvedSCC = ShellResolvedSCC

        self.start("Initialisation_Env")

        # name
        self.set_option("name", image.get_chemical_symbols())
        # number
        self.set_option("number", len(image))

        # type_name
        self.set_option("type_name", set(image.get_chemical_symbols()))
        
        # max_angular_momentum
        maxAngMom = dict()
        for tname in self.type_name:
            maxAngMom.update({tname:get_angular_momentum(tname)})        
        self.set_option("max_angular_momentum", maxAngMom )

        # Pour diminuer le cout de calculs
        if low is False:
            # type_interraction
            interact = lambda atm1,atm2: atm1+"-"+atm2
            self.set_option("type_interraction", 
                np.asarray([ [interact(xj,xi) for xi in self.type_name] 
                                            for xj in self.type_name]).ravel() )
        
        # net_charges
        self.set_option("net_charges",     dict(o=None,l=None,a=None))
        # partial_charges
        self.set_option("partial_charges", dict(o=None,l=None,a=None))

        # angular_momentum
        self._set_angular_momentum()

        # atomic : 
        self._set_atomic_number()
        self._set_atomic_name()
        # azimuthal : 
        self._set_azimuthal_number()
        self._set_azimuthal_name()
        # orbitals : 
        self._set_orbitals_number()
        self._set_orbitals_name()

        # Valences
        self._set_ionic_population()

        self.stop("Initialisation_Env")


    # -------------------------------------------------------------------------
    # Public
    def get_population(self, otype:str="o")->np.ndarray:
        """get_population"""
        return self.population[otype]

    def reshape(self, charges=None):
        count = 0
        resh  = []#np.empty(0, dtype=object)
        for a in self.atomic_number:
            #resh   = np.append(resh, charges[count:count+a] )
            resh.append( charges[count:count+a] )
            count += a
        return np.array(resh,dtype=object) 

    def get_partial_charges(self, otype:str="o")->np.ndarray:
        """get_partial_charges"""
        return self.partial_charges[otype]

    def get_net_charges(self, otype:str="o")->np.ndarray:
        """get_net_charges"""
        return self.net_charges[otype]

    def set_net_charges(self, other:np.ndarray=None, misc:bool=True,)->None:
        """set_net_charges"""

        self.start("Set_Net_Charges")
        n     = np.shape(other)[0]
        table = copy.deepcopy(self.net_charges)
        
        if misc is False:
            table[self.shelltype] = other 
            ptable = copy.deepcopy(table)
            if self.shelltype == "o":
                ptable["o"] -= self.orbitals_valence
            if self.shelltype == "l":
                ptable["l"] -= self.azimuthal_valence
            if self.shelltype == "a":
                ptable["a"] -= self.atomic_valence
        else:
            if n == self.orbitals_name.shape[0]:
                table["o"] = other
                table["l"] = self._reduce( other , otype="l")
                table["a"] = self._reduce( other , otype="a")

            if n == self.azimuthal_name.shape[0]:
                table["o"] = self._extend( other , otype="o")
                table["l"] = other
                table["a"] = self._reduce( other , otype="a")

            if n == self.atomic_name.shape[0]:
                table["o"] = self._extend( other , otype="o")
                table["l"] = self._extend( other , otype="l")
                table["a"] = other

            ptable = copy.deepcopy(table)
            ptable["o"] -= self.orbitals_valence
            ptable["l"] -= self.azimuthal_valence
            ptable["a"] -= self.atomic_valence


        self.set_option("net_charges", table.copy())
        
        self.set_option("partial_charges", ptable)

        #if self.ShellResolvedSCC is False:
        #    self.net_charges     = self.shellResolvedCorrection(self.net_charges)
        #    self.partial_charges = self.shellResolvedCorrection(self.partial_charges)

        #sys.exit(0)
        self.stop("Set_Net_Charges")

    def set_partial_charges(self, other:np.ndarray=None, misc:bool=True)->None:
        """set_net_charges"""

        self.start("Set_Partial_Charges")
        n     = np.shape(other)[0]
        table = copy.deepcopy(self.partial_charges)

        if misc is False:
            table[self.shelltype] = other 
            ptable = copy.deepcopy(table)
            if self.shelltype == "o":
                ptable["o"] += self.orbitals_valence
            if self.shelltype == "l":
                ptable["l"] += self.azimuthal_valence
            if self.shelltype == "a":
                ptable["a"] += self.atomic_valence
        else:
            if n == self.orbitals_name.shape[0]:
                table["o"] = other 
                table["l"] = self._reduce( other , otype="l")
                table["a"] = self._reduce( other , otype="a")

            if n == self.azimuthal_name.shape[0]:
                table["o"] = self._extend( other , otype="o")
                table["l"] = other
                table["a"] = self._reduce( other , otype="a")

            if n == self.atomic_name.shape[0]:
                table["o"] = self._extend( other , otype="o")
                table["l"] = self._extend( other , otype="l")
                table["a"] = other

            ptable = copy.deepcopy(table)
            ptable["o"] += self.orbitals_valence
            ptable["l"] += self.azimuthal_valence
            ptable["a"] += self.atomic_valence

        self.set_option("partial_charges",     table)
        
        self.set_option("net_charges", ptable)

        #if self.ShellResolvedSCC is False:
        #    self.net_charges     = self.shellResolvedCorrection(self.net_charges)
        #    self.partial_charges = self.shellResolvedCorrection(self.partial_charges)


        self.stop("Set_Partial_Charges")

    def shellResolvedCorrection(self, _charge, ):
        """shellResolvedCorrection"""
        #print(_charge)
        for shelltype in ["o","l"]:
            ext = copy.deepcopy(self.angular_momentum)
            if shelltype=="o":
                ext[ext==1.],ext[ext==2.]=4.,9.
            if shelltype=="l":
                ext[ext==2.],ext[ext==1.]=3.,2.

            tmp_atomic,x = np.zeros((int(ext.sum()),))+1.,0
            for i in range(len(ext)):
                t = int(ext[i])
                tmp_atomic[x:x+t] *= _charge["a"][i] / t
                x += t
            _charge[shelltype] = tmp_atomic.copy()
        #print(_charge)
        return _charge


    # -------------------------------------------------------------------------
    # Private

    def _set_atomic_number(self,):
        """_set_atomic_number"""
        mam   = self.angular_momentum
        other = np.array([_get_momentum(mom) for mom in mam],dtype=np.int32)
        self.set_option( "atomic_number",other)    

    def _set_azimuthal_number(self,):
        """_set_azimuthal_number"""
        other = []
        for i,elm in enumerate(self.atomic_number):
            if elm>=1:_tab=[1]
            if elm>=4:_tab+=[3]
            if elm>=9:_tab+=[5]
            other += _tab
        self.set_option( "azimuthal_number",np.array(other,dtype=np.float64))

    def _set_orbitals_number(self,):
        """_set_orbitals_number"""
        other = []
        for i,elm in enumerate(self.atomic_number):
            if elm>=1:_tab=[1]
            if elm>=4:_tab+=[3]*3
            if elm>=9:_tab+=[5]*5
            other += _tab
        self.set_option( "orbitals_number",other)    

    def _set_atomic_name(self,):
        """_set_atomic_name"""
        other = np.arange(0, self.number )
        self.set_option( "atomic_name",other)    

    def _set_azimuthal_name(self,):
        """_set_azimuthal_name"""
        count=-1
        other=np.zeros_like(self.azimuthal_number)
        for i,elm in enumerate(self.azimuthal_number):
            if elm == 1:count+=1;other[i]=count
            if elm == 3:other[i]=count
            if elm == 5:other[i]=count
        self.set_option( "azimuthal_name",other)    

    def _set_orbitals_name(self,):
        """_set_orbitals_name"""
        count=-1
        other=np.zeros_like(self.orbitals_number)
        for i,elm in enumerate(self.orbitals_number):
            if elm == 1:count+=1;other[i]=count
            if elm == 3:other[i]=count
            if elm == 5:other[i]=count
        self.set_option( "orbitals_name",other)        
        


    def _set_ionic_population(self):
        """_set_ionic_population"""

        self.start("Set_Ionic_Charges",)
        _valence = []
        for i,sym in enumerate(self.image.get_chemical_symbols()):
            _valence.append( list(atomic_data[sym].atom_studied["occupations"].values()) )
        mam = self.angular_momentum

        atm_valence,orb_valence = [0.]*len(_valence),[]
        for i,elm in enumerate(_valence):
            atm_valence[i],tmp = np.sum(elm),[]
            for j,val in enumerate(elm):
                if len(tmp) >= mam[i]+1:continue
                if j==0:tmp.append(val)
                if j==1:tmp.append(val)
                if j==2:tmp.append(val)
            orb_valence+=tmp
        
        atm_valence=np.array(atm_valence,dtype=np.float64).ravel()
        self.set_option("atomic_valence",atm_valence)
        
        orb_valence=np.array(orb_valence,dtype=np.float64).ravel()
        self.set_option("azimuthal_valence",orb_valence)

        self.set_option("orbitals_valence",self._extend(orb_valence,otype="o"))
    
        self.stop("Set_Ionic_Charges")

    def _set_distances_nopbc(self, units=1.):
        pos = self.image.positions * units
        Distance_matrix = scipy.spatial.distance.cdist(
                pos, pos, 'euclidean')
        self.set_option("distance_matrix",Distance_matrix)

    def _set_interraction_nopbc(self):
        """_set_interraction_nopbc"""
        type_name = list(self.name)
        interact = lambda atm1,atm2: atm1+"-"+atm2
        self.set_option("interraction_matrix", 
            np.asarray([ [interact(xj,xi) for xi in type_name] 
                                          for xj in type_name]) )
    
    def _set_angular_momentum(self):
        """_set_angular_momentum"""
        element = self.name
        angM = np.array([get_angular_momentum(elm) for elm in element],dtype=np.float64)
        self.set_option( "angular_momentum", angM)







    # -------------------------------------------------------------------------
    # Utils

    def _reduce(self, other, otype="l"):
        """_reduce"""
        #assert len(other)>len(self.atomic_number), TypeError()
        tmp  = 0
        save = [] 
        if otype == "l":
            number = self.azimuthal_number
            name   = self.azimuthal_name
        elif otype == "a":
            number = self.atomic_number
            name   = self.atomic_name
        else:
            raise TypeError("Reduction of charge isn't allowed.")

        for i,(num,nam) in enumerate(zip(number,name)):
            tmp   += num
            lp,lm  = int(tmp),int(tmp-num)
            save.append(other[lm:lp].sum())
        return np.array(save,dtype=np.float64)

    def _extend(self, other, otype="l"):
        """_extend"""
        save = []

        if len(other)==len(self.atomic_name):
            if otype=="l":l=[1,2,3]
            if otype=="o":l=[1,4,9]
            number = self.atomic_number
            name   = self.atomic_name
            for i,(num,nam) in enumerate(zip(number,name)):
                if num == 1:save+=[other[i]]
                if num == 4:save+=[other[i]/l[1]]*l[1]
                if num == 9:save+=[other[i]/l[2]]*l[2]

            return np.array(save,dtype=np.float64)

        if len(other)==len(self.azimuthal_name):
            number = self.azimuthal_number
            name   = self.azimuthal_name
            for i,(num,nam) in enumerate(zip(number,name)):
                if num == 1:save+=[other[i]]
                if num == 3:save+=[other[i]/num]*int(num)
                if num == 5:save+=[other[i]/num]*int(num)

            return np.array(save,dtype=np.float64)

