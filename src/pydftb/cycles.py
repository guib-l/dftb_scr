#!/usr/bin/env python3
import __future__

import os
import sys
import copy
import time

import scipy
import numpy as np

import ase
from ase import Atoms


class AbstractCycle(object):
    def __init__(self, verbose = True, filename="dftb.output"):
        """ __init__ """
        self._attach_function = []
        self.verbose = verbose

        self.filename = filename

        if os.path.exists(filename):
            print(" [WARNING] Erase file name  : ",filename,'\n')
            os.remove(filename)  

    def _generate_line(self, *args):
        """
        _generate_line

        Fonction qui sert à générer une chaine de caratère (string) à partir 
        d'un tableau dans un format scientifique.
        """
        line = ""
        for q in args:
            if isinstance(q, str): line += "{}".format(q)
            elif isinstance(q, int): line += "{}\t".format(q)
            else: line += "{:1.8e}\t".format(q)
        return line

    def _write_convergeance(self, *args):
        """
        _write_convergeance
        
        Fonction qui écrite dans un fichier une chaine de caractère générée avec 
        '_generate_line'. Autant d'éléments peut lui être passé.
        """
        line = ""
        with open(self.filename,"a") as fin :
            line = self._generate_line(*args)
            line += "\n"
            fin.write(line)

    def attach(self, function, period=1, *args):
        """attach"""
        assert callable(function), TypeError
        self._attach_function.append( {"function":function,"args":args,"cycle":period,} )

    def iteration(self,):
        """iteration"""
        raise NotImplementedError

    def run(self,):
        """run"""
        raise NotImplementedError














