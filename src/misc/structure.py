#!/usr/bin/env python3
import os
import sys

import numpy as np

import ase
from ase import Atoms
import typing

from .plot import draw_dictionnary

class structure(object):
    def __init__( self, 
                  atoms:ase.Atoms, 
                  **option:dict):
        """
        :: structure ::

        Classe qui permet de mettre en œuvre une structure de donnée claire et 
        inter-opérable dans tout le code. 
        Ce morceau de code ne peut renvoyer un résultat qu'à la condition que 
        le dictionnaire 'results' soit remplie d'élément non nuls. Cela à pour 
        de contraindre le système planter e cas de mauvaise entrée plutôt que
        de continuer et faire n'importe quoi.

        args :   
            atoms (ase.Atoms) : 
                Atoms object de la librairie ase.

            options (dict) : 
                dictionnaire qui sera entré dans la structure en plus des élements 
                standard tels que l'atome ou le poids.

        return :
            atoms (ase.Atoms) :
                Retourne un Atoms object de la librairie ase.

            weight (float) :
                Retourne un flottant défini par défaut à 1.

            get_results (dict) : 
                Retourne un dictionnaire qui va être composé des éléments 'atom',
            ainsi que 'weight' en plus des éléments qui lui auaront été ajouté par 
            l'utilisateur.

            reset (None) :
                Vide tout les éléments sauf 'atom' et 'weight' et les fixe à None.

            set_empty (dict) :
                Retourne un dictionnaire vidé de son contenu sauf des clefs.

        use :
        
        >>> from pylib.series import structure
        >>> A = structure(atoms , **data)
        
        """
        assert isinstance(atoms, ase.Atoms), \
                        ValueError("atoms needs to be a ase.Atoms object!")
        self._results = {"atom":atoms, "weight":1.}

        for key in option.keys():
            self.set_option(key, option[key])
            
    @property
    def atoms(self)->ase.Atoms:
        return self._results["atom"]

    @property
    def weight(self)->float:
        return self._results["weight"]

    def __check(self)->bool:
        """ Méthode qui check si les élément diu dictionnaire sont nuls ou non. """
        state = True
        for key in self._results.keys():
            if self._results[key] is None:
                state=False
        return state

    def _check_key(self,key:str)->bool:
        state = False
        for k in self._results.keys():
            if key.lower() == k.lower():
                state = True
        return state

    def set_empty(self)->dict:
        """ Retourne un dictionnaire vidé de son contenu sauf des clefs. """
        tmp=dict()
        for key in self._results.keys():
            tmp[key] = {}
        return tmp   

    def reset(self)->None:
        """ Renvoie tout les éléments sauf 'atom' et 'weight' à None. """
        for key in self._results.keys():
            if key=="atom":continue
            if key=="weight":continue
            self._results[key] = None

    def set_option( self, 
                    name:str, 
                    obj)->None:
        """
        Initialise un nouvel élément au dictionnaire de manière à ce qu'il 
        non nul. 
        """
        self._results.update({name:obj})
        if self._check_key(key=name) is False:
            raise ValueError("Some elements { %s } are set to 'NoneType'")
        
        setattr(self, name, obj)

    def get_results(self, alert=False )->dict:
        """
        Retourne le résultat.
        """
        if alert is not False:
            if self.__check() is False:
                print("[ALERT] :: Some elements are set to 'NoneType'")
        return self._results

    def __call__(self, key):
        return self._results[key]

    def __repr__(self):
        return draw_dictionnary(self._results)

            