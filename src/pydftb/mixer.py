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

from atomic.atomic import Atomic
from atomic.environnement import Environnement



class MixCycle(object):
    def __init__(self, method="simple", **history):
        """__init__"""
        self.history = history
        self.history.update(saved = [0,] )

        if method == "simple":
            self._mixer_call = self._simple

    def _simple(self, x):
        """_simple"""
        x0               = self.history['x0']
        mixing_parameter = self.history['mixing_parameter']

        # Sauvegarde x0
        self.history["saved"].append(x0)
        
        _news = x * mixing_parameter + x0 * (1 - mixing_parameter)
        
        # Update x0
        self.history.update(x0 = _news)
        return _news

    def mixing(self, x=None, **kwargs):
        """mixing"""
        return self._mixer_call(x, **kwargs)




