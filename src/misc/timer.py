import os
import sys
import math
import time
import random
import numpy as np

import matplotlib.pyplot as plt


class Timer:
    def __init__( self, label:str)->None:

        self.label    = label
        self.timer_id = {}
        self._start   = {}
        self._stop    = {}
        self._level   = {}

        tmp = {label:time.time()}
        self._start.update(tmp)
        self.is_running = 1

    def reset_timing(self,):
        self.timer_id = {}
        self._start   = {}
        self._stop    = {}
        self._level   = {}

        tmp = {self.label:time.time()}
        self._start.update(tmp)
        self.is_running = 1

    def sumerize(self):
        self.is_running = 0
        l=0
        print("\n : TIMER : ",self.label)
        for timer in  self.timer_id.keys():
            old_l = l
            l = self._level[timer]
            print("{}- {:15} = {}".format("\t"*l,timer,float(self.timer_id[timer])))
            print("\t---------------")
        print(" : ")
        return self.timer_id

    def start(self,new_instance, level=1):
        if self.is_running == 1:
            tmp = {new_instance:time.time()}
            self._start.update(tmp)
            tmp = {new_instance:level}
            self._level.update(tmp)
        
    
    def stop(self,instance="instance"):
        if instance in self._start:
            tmp = {instance:time.time()}
            self._stop.update(tmp)
            tmp = {instance:self._stop[instance] - self._start[instance]}
            self.timer_id.update(tmp)
        else:
            raise AssertionError('Timer %s cannot be stopped; it is not running!' % self.label)

        if instance == self.label:
            print("Timer {} close;".format(self.label))
            self.is_running = 0


