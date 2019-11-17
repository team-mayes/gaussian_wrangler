import jpype
import jpype.imports
from jpype.types import *




class HartreeWrapper():
    def __init__(self):
        jpype.addClassPath("hartree/*")
        jpype.startJVM(convertStrings=False)

        from org.cmayes.hartree.loader.gaussian import SnapshotLoader
        from java.lang import System
        from java.io import FileReader

        self.loader = SnapshotLoader()

    def read_gaussian(self, files):
        mapped_results = {}
        for cur_file in files:
            mapped_results[cur_file] = self.loader.load(cur_file, FileReader(cur_file))

        return mapped_results

    def __del__(self):
        jpype.shutdownJVM()


