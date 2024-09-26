#!/usr/bin/env python3
import os
import sys

import numpy as np

import ase
from ase import Atoms
import typing

from .plot import draw_dictionnary
from .structure import structure

        
class series(object):
    def __init__(self, 
                 name:str="structure_series_00",
                 series:typing.Union[typing.Dict, typing.List ]=None)->None:
        """
        :: series ::

        Classe qui permet de mettre en œuvre une structure de donnée claire et 
        inter-opérable dans tout le code. Object qui se compose comme une liste 
        de plusieurs structures. Ce dernier peut être employé comme une liste.

        args :   
            name (str) :
                Nom de la series de structure.
            series (structure) :
                Tableau de structure.

        return :
        """
        self._series = []
        self.idx = 0
        if series is not None:
            self._series = series

    def __call__(self)->typing.Union[typing.Dict, typing.List ]:
        """ Retourne la liste de structure """
        return self._series

    def __len__(self)->int:
        """ Longeur de la liste """
        return len(self._series)

    def __getitem__(self, index:int)->structure:
        """ Retourne un élement de la liste """
        return self._series[index]

    def __iter__(self):
        """ Itération """
        for i in self._series :
            yield i

    def __next__(self)->structure:
        """ Itérateur """
        try:
            item = self._series[self.idx]
        except IndexError:
            raise StopIteration()
        self.idx += 1
        return item

    def __del__(self)->None:
        """ Vide la liste """
        self._series = []
        self.idx = 0

    def __new__(cls, *args, **kwargs):
        instance = super(series, cls).__new__(cls, *args, **kwargs)
        return instance

    def __add__(self, obj:structure)->None:
        """ Ajoute une structrue """
        assert isinstance(obj, structure), "Only accept structure object !"
        self._series.append(obj)

    def __check(self, props:str)->bool:
        """ Vérifie l'existence d'une propriété """
        for i,struc in enumerate(self._series):
            try:
                tmp = struct[props]
                state = True
            except:
                state = False
                break
        return state

    def add(self, obj:structure)->None:
        """ Ajoute une structrue """
        assert isinstance(obj, structure), "Only accept structure object !"
        self._series.append(obj)

    def emptySeries(self)->typing.Union[typing.Dict, typing.List ]:
        """ 
        TODO : Peut être amélioré.  (data.keys())
        """
        return series( name="< Empty.Series.Name >",
                       _series=[{i["atom"],i["weight"],} for i in self._series])

    def copy(self)->typing.Union[typing.Dict, typing.List ]:
        """
        """
        return self.emptySeries()
    
    def __repr__(self):
        return "".join(self._series.__str__()) 
    
    def __str__(self):
        return "".join(self._series.__str__()) 

    