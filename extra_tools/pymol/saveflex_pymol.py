from pymol import stored
from pymol import cmd
import numpy as numpy


def saveflex_from_mol2 ( Filename, userSelection ):
# saves bonds within selection in .flex format

        # get atoms indexes / IDs from pymol selection

        # IF PROBLEMS CHECK IF MOL2 WAS GENERATED USING INDEX AS KEY
        # index = numpy.array( cmd.index( userSelection ) )    
        # index_list = [ int(index[k, 1]) for k in range( len(index) ) ]


        index_list = numpy.array( cmd.identify( userSelection ) )

        # open needed files
        molFilename = Filename + ".mol2"
        flexFilename = Filename + ".flex"
        mol2 = open(molFilename, "r")
        flex = open(flexFilename, "w")

        print(index_list)

        # gets bond list from mol2 and prints selected flexible bonds ( index - 1 )
        bondsection = 0
        line = mol2.readlines()
        for i in line:

                currentLine = i.split()

                if len(currentLine) > 0:

                        if currentLine[0] == "@<TRIPOS>BOND" :
                                bondsection = 1
                        elif currentLine[0] == "@<TRIPOS>SUBSTRUCTURE" :
                                bondsection = 0
                                break

                        elif bondsection:
                            if int( currentLine[1] ) in index_list and int( currentLine[2] ) in index_list:
                                    print >> flex, int(currentLine[1]) - 1, int(currentLine[2]) - 1
                                    # prints in pymol console
                                    print(int(currentLine[1]) - 1, int(currentLine[2]) - 1)

                                                                                                                         
