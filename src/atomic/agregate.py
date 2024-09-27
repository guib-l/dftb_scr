#!/usr/bin/env python
import os
import io
import sys
import math
import time
import copy
import random
import numpy as np
from copy import deepcopy

import ase
from ase import Atoms
from ase.optimize import QuasiNewton,LBFGS
from ase.calculators.dftb import Dftb

from sklearn.metrics.pairwise import euclidean_distances

from src.misc.storage import Storage

class Agregate(object):
    def __init__( self, 
                  atoms          = None,
                  tip_name       = "Tip-00",
                  tip_center     = [0.,0.,0.], 
                  tip_dimension  = { "R":5., }, 
                  clip_agregat         = True,
                  monteCarlo_atoms     = False,
                  monteCarlo_positions = False,
                  monteCarlo_params = {
                      "number"     : 0,
                      "noise"      : 0.00,
                       "centering" : True},
                  seed=None,
                  *args, 
                  **kwargs):

        super().__init__(*args, **kwargs)

        if isinstance(atoms, str):
            tip_name = atoms.split('.pkl')[0]
            
            atoms = Storage(atoms).stream()

            assert isinstance(atoms, ase.Atoms), \
                    TypeError("atoms isn't a ase.Atoms object yet.")
            tip_center    = atoms.get_center_of_mass()
            tip_dimension = { "R":1e6, }

        self._tip_center    = tip_center
        self._tip_dimension = tip_dimension


        center_of_masses = atoms.get_center_of_mass() - self._tip_center
        atoms.positions -= center_of_masses
        atoms.cell = None

        self.tip_name = tip_name
        
        self.seed = None
        if seed is not None:
            self.seed = seed

        if clip_agregat:

            # Choix primaire du type de sonde
            if len( list(self._tip_dimension.keys()) ) == 2:
                
                self.atoms = self._cylindre(atoms)

            elif np.all( list(self._tip_dimension.keys()) == np.asarray(["X","Y","Z"]) ):

                self.atoms = self._squared(atoms)
            else:
                self.atoms = self._rounded(atoms)

        else:
            self.atoms = atoms.copy()

        if monteCarlo_atoms:
            self.atoms = self._random_remove_atoms(
                **monteCarlo_params)
        if monteCarlo_positions:
            self.atoms = self._random_noising_positions(
                **monteCarlo_params)
        else: pass

    
    def _random_noising_positions(self, noise=0.01, centering=True, **kwargs):
        
        atoms = self.atoms.copy()
        if self.seed is not None:
            rng = np.random.default_rng(seed=self.seed)
        else: rng = np.random.default_rng()

        if not centering:
            y_noise = np.asarray(
                [[noise, noise, noise] for x in atoms.positions]
            ).ravel() * rng.normal(size=atoms.positions.size)

        if centering:

            l = np.sqrt(np.sum(
                (atoms.positions - atoms.get_center_of_mass())**2,
                axis=1))
            fact = (l / np.max(np.abs(l)))**2
            y_noise = np.asarray(
                [[noise, noise, noise] for x in atoms.positions]
            ).ravel() * rng.normal(size=atoms.positions.size)

        
        noise_r = np.reshape(
            y_noise,(len(atoms),3)
        ) * fact.reshape((-1,1))
        atoms.positions += noise_r
        return atoms.copy()



    def _random_remove_atoms(self, number, **kwargs):
        """_random_remove_atoms"""
        assert number < len(self.atoms)

        atoms = self.atoms.copy()
        if self.seed is not None:
            rng = np.random.default_rng(seed=self.seed)
        else: rng = np.random.default_rng()

        for i in range(number):

            numR = int(rng.uniform(0, len(atoms)))
            atoms.pop(numR)

        return atoms.copy()


    def get_tip_state(self):
        return dict(
            tip_name      = self.tip_name,
            tip_center    = self._tip_center,
            tip_dimension = self._tip_dimension
        )

    def _squared(self, atoms):
        """_squared"""
        symbols = atoms.get_chemical_symbols()

        new_symbols   = []
        new_positions = []

        for (sym,pos) in zip(symbols,atoms.positions):

            if pos[0] > -self._tip_dimension["X"]/2 and \
               pos[0] < self._tip_dimension["X"]/2 :
                
                if pos[1] > -self._tip_dimension["Y"]/2 and \
                   pos[1] < self._tip_dimension["Y"]/2 :
                    
                    if pos[2] > -self._tip_dimension["Z"]/2 and \
                       pos[2] < self._tip_dimension["Z"]/2 :
                        
                        new_positions.append(pos)
                        new_symbols.append(sym)

        atoms = Atoms( "".join(new_symbols), 
                       positions=new_positions,
                       cell=None,
                       pbc=False) 

        return atoms


    # Calcul une sphère autour du centre de masse de la super-cellule
    def _rounded(self, atoms):
        """
        _rounded :: 

            Fonction qui définie une sphère de rayon R. CE sera la seule chose 
        qui lui sera absolument nécessaire. Pour l'instant, on ne peut pas 
        wrapper les structures entre elles. On spécifie une sphère d'atome de 
        la sorte :
            { "R" : (float) }

        args : 
            atoms (ase.Atoms) :
            
            Object ase.Atoms avce les types chimique et les positions atomique 
        absolues. 

        return :
            atoms (ase.Atoms) : 

            Nouvel object ase.Atoms avce les types chimique et les positions 
        atomique absolues sous la forme d'une sphère. 

        """
        symbols = atoms.get_chemical_symbols()

        new_symbols   = []
        new_positions = []

        for (sym,pos) in zip(symbols,atoms.positions):
        
            if np.sqrt(np.sum(pos**2)) <= self._tip_dimension["R"]:
                new_positions.append(pos)
                new_symbols.append(sym)

        atoms = Atoms( "".join(new_symbols), 
                       positions=new_positions,
                       cell=None,
                       pbc=False) 

        return atoms

    def _cylindre(self, atoms):
        """
        _cylindre :: 

            Fonction qui définie un cylindre de rayon R dans l'axe choisi. 
        Ce sera la seule chose qui lui sera absolument nécessaire (un rayon 
        et un axe). Pour l'instant, on ne peut pas wrapper les structures 
        entre elles. On spécifie une sphère d'atome de 
        la sorte :
            {   "R" : (float),
                "X" : (float)  -> A choisir entre (X,Y,Z) }

        args : 
            atoms (ase.Atoms) :
            
            Object ase.Atoms avce les types chimique et les positions atomique 
        absolues. 

        return :
            atoms (ase.Atoms) : 

            Nouvel object ase.Atoms avce les types chimique et les positions 
        atomique absolues sous la forme d'un cylindre. 

        """
        symbols = atoms.get_chemical_symbols()

        new_symbols   = []
        new_positions = []

        ref = ["X","Y","Z"]
        new_axes = []
        for i,axe in enumerate(ref):            
            if axe not in self._tip_dimension: new_axes.append(i) 
            else : atribut = (i,axe)

        assert len(new_axes)>=2, ValueError()

        for (sym,pos) in zip(symbols,atoms.positions):
    
            npos = np.asarray ( [pos[ new_axes[0]], pos[new_axes[1] ] ])

            if np.sqrt(np.sum(npos**2)) <= self._tip_dimension["R"]:

                if (pos[atribut[0]] < self._tip_dimension[atribut[1]]/2) and \
                   (pos[atribut[0]] > -self._tip_dimension[atribut[1]]/2):

                    new_positions.append(pos)
                    new_symbols.append(sym)


        atoms = Atoms( "".join(new_symbols), 
                       positions=new_positions,
                       cell=None,
                       pbc=False) 

        return atoms






