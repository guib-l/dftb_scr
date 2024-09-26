#!/usr/bin/env python3
import __future__

# Import standard de python3
import os
import sys
import numpy as np

STEP_verbose_debug    = 1000
STEP_verbose_path     = 100
STEP_verbose_standard = 10
STEP_verbose_low      = 1
STEP_verbose_not      = 0

# *************************** \
def progressbar(it, prefix="", size=80, out=sys.stdout): # Python3.3+
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print("{}[{}{}] {}/{}".format(prefix, u'█'*x, "."*(size-x), j, count), 
                end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)











