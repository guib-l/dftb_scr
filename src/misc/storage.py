import sys
import os
import pickle


class Storage:
    def __init__( self, 
                  filename:str, 
                  txt:_io.TextIOWrapper=sys.stdout)->None:
        """
        ::  IOmodel ::

        Classe qui peut lire et enregistrer 
        des objets dans des fichiers en binaire avec l'emploie de 
        la librairie Pickle.
        
        """
        self._txt = txt
        
        self._filename = filename
        self._model = None

    def read(self)->None:
        """ Méthode pour lire les fichier """
        self._model = pickle.load(open(self._filename,'rb'))        

    def write( self, 
               other, 
               anot:str='wb')->None:
        """ Méthode pour écrire dans un fichier """

        pickle.dump(other, open(self._filename, anot) )

    def stream(self):
        """ Méthode qui retourne depuis un fichier un objet 'model'. """
        if self._model is None:
            self.read()
        return self._model
        
