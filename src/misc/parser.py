import os
import sys

import numpy as np

class readTable :
    def __init__( self, 
                  file_name:str,
                  comment:str = "",
                  verbose:bool = False,
                  replace:float=None,):
        

        self.file_name  = file_name
        self.comment    = comment
        self._verbose   = verbose
        self.data       = None
        
        self.__call__(replace)

    def __collect(self, replace=None):
        fichier = open(self.file_name, "r")
        assert fichier
        
        lst = []
        for lines in fichier:
            if lines[0] != "#" :
                tmp = lines.rstrip('\n\r').split()
                for i in range(len(tmp)):
                    if tmp[i]=='':pass
                    else :
                        try:tmp[i]=float(tmp[i])
                        except:pass
                lst.append(tmp)
        fichier.close()

        row_lengths=[]
        for row in lst:
            row_lengths.append(len(row))
        max_length = max(row_lengths)

        for row in lst:
            while len(row) < max_length:
                row.append( replace )

        self.data = lst
        
    def print(self):
    	if self.verbose == True:
    	    print("debug :: Les donnée sont issues du fichier {}"
                  .format(self.file_name))
    	    print("debug :: Commentaire :: {}"
                  .format(self.comment))
    	else:
            print("Need to allow verbose variable")
            
    def extract(self, column_1):
    	return np.array(self.data[:, column_1])

    def __call__(self, replace=None):
    	self.__collect(replace=replace)
    	if self._verbose == True:
            self.print()

    def write( self, 
               newFile:str = "", 
               comment:str = ""):

        M = self.data.shape[1]
        N = len(self.data[0])
        
        fw = open(newFile,"w")
        fw.write("# " + comment + "\n")
        for j in range(N):
            for i in range(M):
                fw.write("{}".format(self.data[i][j]))
                if i < M-1:fw.write("\t")
            fw.write("\n")
        return 0




