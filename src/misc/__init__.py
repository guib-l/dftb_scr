#!/usr/bin/env python3
import __future__

# Import standard de pythons3
import os
import sys
import re
import math
import numpy
import scipy


from data.misc import progressbar
from data.parser import readTable
from data.picklemodel import pickelFormat,Storage
from data.plot import ( draw_dictionnary, Plot )
from data.series import ( series, structure )
from data.timer import Timer

