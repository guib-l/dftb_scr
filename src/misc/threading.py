#!/usr/bin/env python
import os
import io
import sys
import math
import time
import random
import numpy as np

import queue
import threading


STEP_verbose_debug    = 1000
STEP_verbose_path     = 100
STEP_verbose_standard = 10
STEP_verbose_low      = 1
STEP_verbose_not      = 0


# Pas sûr de ce morceau de code


class IsolatedThread (threading.Thread):

    def __init__(
            self, 
            threadID, 
            threadName, 
            process,
            argument,
            resultsList,
            txt=sys.stdout,
            verbose=0):
        
        threading.Thread.__init__(self)

        self._threadID = threadID
        self._name = threadName
        self._queue = process
        self._args = argument
        self._results = resultsList

        self.txt = txt
        self.verbose = verbose
        
    def run(self):
        if self.verbose > STEP_verbose_standard:
            print(" > [WARNING] Starting " + self._name, file=self.txt)

        self.process_data()

        if self.verbose > STEP_verbose_standard:
            print(" > [WARNING] Exiting " + self._name, file=self.txt)

    def process_data(self):
        
        while not exitFlag:
            queueLock.acquire()

            if not workQueue.empty():
                data = self._queue.get()
                arg  = self._args.get()
                queueLock.release()

                if self.verbose > STEP_verbose_standard:
                    print(" > [ALERT] %s processing\n\t %s\n" % (self._name, data))

                try : 
                    self._results.append( data(arg) ) 
                    #time.sleep(0.1)
                
                except : pass

            else:
                queueLock.release()


class ConcatenateThread(object):

    def __init__(
            self,
            nThreads,
            processusList,
            argumentList,
            txt=sys.stdout,
            verbose=100,
            **kwargs):

        assert isinstance(nThreads, int), \
                ValueError(" > [ERROR] Number of threads is not an 'int' !")
        self._nThreads = nThreads
        
        self.txt = txt
        self.verbose = verbose
        
        self._processus = processusList
        self._argument  = argumentList
        self._workQueue = queue.Queue(0)
        self._argsQueue = queue.Queue(0)
        self._output = []

        self.__fill_queue()
        self.__start_threads()

        
    def __fill_queue(self):

        # Fill the queue
        queueLock.acquire()

        for word in self._processus:
            workQueue.put(word)

        for arg in self._argument:
            argsQueue.put(arg)

        queueLock.release()

        self._queueLock = queueLock

    def __start_threads(self):

        self.threads = []

        # Create new threads
        for iDThread in progressbar(range(self._nThreads), "Threads init  : ", 60):
            #for iDThread in range(self._nThreads):
            
            tName = "thread-%s_%s"%(iDThread+1,time.time())

            thread = IsolatedThread(
                iDThread+1, 
                tName, 
                workQueue, 
                argsQueue,
                self._output )

            thread.start()
            self.threads.append(thread)
                
        # Wait for queue to empty
        while not workQueue.empty():
            pass


    def run_threads(self):

        # Wait for all threads to complete
        #for t in self.threads:
        for t in progressbar(range(len(self.threads)), "Threads prog  : ", 60):

            self.threads[t].join(timeout=None)

        #print("Exiting Main Thread")


class MultiThread(ConcatenateThread):
   

    def __init__(
            self,
            nThreads,
            processusList, 
            argument,
            **kwargs):
            
        global exitFlag
        global workQueue,queueLock, argsQueue

        # Code d'exécution de sortie des threads
        exitFlag = 0
        # Génère une queue infinie
        workQueue = queue.Queue(0)
        # Génère une queue infinie
        argsQueue = queue.Queue(0)
        # Défini un veroux sur chaque thread de sorte qu'ils soient synchrones.
        queueLock = threading.RLock()

        super().__init__(
                nThreads,
                processusList, 
                argument, 
                **kwargs)

        exitFlag = 1

        self.run_threads()

    def get_output(self):

        return self._output


if __name__ == '__main__':

    fct3 = lambda x: x**4
    fct2 = lambda x: x**3
    fct = lambda x: x**2
    nameList = [fct, fct, fct2,fct2,fct3, fct, fct, fct2,fct2,fct3]

    #exitFlag = 0
    M = MultiThread(2, nameList)

    #exitFlag = 1

    #M.run_threads()

    print(M._output)
