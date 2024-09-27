#!/usr/bin/env python3
import __future__

# Import standard de pythons3
import os
import sys
import copy
import time

import scipy
import numpy as np

import ase
from ase import Atoms

# ******************************** #

from atomic.atomic import Atomic

from .mixer import MixCycle
from .cycles import AbstractCycle


class DftbSCC(AbstractCycle, Atomic, MixCycle):
    """SCCcycle"""
    def __init__(self, 
            image:ase.Atoms  = None,
            method           = "simple", 
            mixer_argument   = {
                "mixing_parameter" : 0.1,},
            verbose = True,
            **kwargs )->None:
        """ __init__ """
        AbstractCycle.__init__(self, verbose = verbose)
        Atomic.__init__(self, image, **kwargs)
        MixCycle.__init__(self, method=method, **mixer_argument)

        self.q0, self.q1 = None, None

    def iteration(self, charges=None, temperature=None):
        """iteration"""
        self.diagonalize(charges=charges, eigMethod=scipy.linalg.eigh)
        self.set_population(temperature=temperature)
        self.set_charges()

    def callback(self,):
        print( " Energie total : ", self.energies )
        print( " ............... ", self.get_total_energy())       
        print( " ............... ", self.diff_energy)      


    def run(self,
            epoch:np.int32         = 3,
            temperature:np.float64 = 100.,
            tolerance:np.float64   = 1e-8,
            **kwargs )->None:
        """run"""
        if self.verbose:
            print("\n","-"*80)
            print(" > Start calculation  DFTB : ")
            print("-"*80,"\n")
        self.step = 0

        old_energy = 0.
        new_energy = 0.

        self.q0            = self.get_partial_charges()
        self.history['x0'] = self.q0

        while True: 
            if self.verbose:
                print("\n Iteration     :  %s"%self.step)

            self.iteration( self.q0, temperature=temperature )

            self.q1 = self.get_partial_charges()
            self.set_total_energy(charges=self.q1)

            new_energy = self.get_total_energy()
            self.diff_energy =  new_energy - old_energy 
            old_energy = new_energy

            for elm in self._attach_function:
                if self.step % elm["cycle"] == 0:
                    elm["function"]( *elm["args"] )

            self.q0 = self.mixing(self.q1)

            self.set_partial_charges(self.q0)
            self.reset()

            if self.step >= epoch: break
            if np.abs(self.diff_energy) <= tolerance: break
            
            self.step += 1

        return 0
