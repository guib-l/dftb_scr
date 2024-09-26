#!/usr/bin/env python3
import __future__

# Import standard de python3
import os
import sys
import numpy as np


class Plot(object):

    def __init__( 
            self, 
            *args, 
            indentation=0, 
            separateur="-", 
            **kwargs ):

        self.__indentation = indentation
        self.__separateur = separateur

        super().__init__(*args, **kwargs)

    def plot_prototype(self):
        cls = type(self)
        txt = None#sys.stdout
        separateur = self.__separateur * 80
        
        print(separateur, file=txt)
        print("Class Object name :", cls.__name__, file=txt)
        print(" > DOC : ", cls.__doc__, file=txt)
        print(" > MRO : ", cls.mro(), file=txt)

        print(" > Arguments : ", file=txt)

        indentation = ' ' * self.__indentation
        
        for attr in self.__dict__:
            if not attr.startswith('_'):
                valeur = getattr(self, attr)
                print(f"\t{indentation}{attr} = {valeur}", file=txt)

        print(separateur, file=txt)


# *************************** \
def __print_dict( other:dict, 
                  space:int=0,
                  indent:str="\n",
                  txt=sys.stdout):
    results = "{}".format(indent)
    for keys in other.keys():
        space+=1
        if isinstance(other[keys],dict):
            results += "{} * {}  {}".format(" "*space, keys,indent)
            __print_dict(other[keys],space,txt)
        else:
            results += "{}{} {}{}".format("  "*(space-1)+"| ", keys,other[keys],indent)
        
        space-=1
    return results

# *************************** \
def draw_dictionnary( other:dict,
                      name="Print dictionnary",
                      indent="\n",
                      txt=sys.stdout):
    space=0
    result = __print_dict(other, space, indent, None)
    return "{} *** {}\n{}\n *** {}".format(indent,name,result,indent)
