# Solve a mixed time-energy optimal control problem in 2 dimensions and
# saves the result

# Nicholas Sharp - nsharp3@vt.edu

import os,sys,time
from optparse import OptionParser
import shutil

import numpy as np

from AlphaHullProp_2D import AlphaFrontPropTESolver
from FluidFuncs import *


def main(argV):

    ## Parse options
    usage = '''python CurrentControlSolver_2D.py INPUT_ARGS_FILE \n\n \
Run the 2D time-energy optimal alpha shape solver'''

    parser = OptionParser(usage=usage)
    parser.add_option("", "-o", help="Output directory", metavar="string", default='./')

    (opts, args) = parser.parse_args(argV)

    if len(args) < 1:
        print("Parameters file required")
        exit()

    print('Using output directory: "%s"'%(opts.o))

    ## Read parameters from file
    '''
    Type:\tTYPE
    Flow Name:\tNAME
    Flow Params:\tp1\tp2\t...
    p0:\tp0
    pf:\tpf
    rhoE:\trhoE
    rhoT:\trhoT
    delJ:\tdelJ
    Initial Speed Upper Bound:\tspeed
    Completion Distance:\tdist
    Alpha:\talpha
    '''
    fileName = args[1]
    print("Reading parameters from file: " + fileName)
    # Get all lines without comments and copy the output file
    lines = []
    for line in open(fileName,'r').readlines():
        if line[0] is not '#':
            lines.append(line.strip())

    # Save a copy of the parameter file
    shutil.copyfile(fileName, opts.o + "input_params.txt")

    # Verify that this is the correct run type
    runType = lines[0].split("\t")[1]
    if runType != "2D-TIME-ENERGY":
        print("Incorrect run type, must be 2D-TIME-ENERGY")
        exit()

    # Get the type of flow and flow parameters
    flowName = lines[1].split("\t")[1]
    flows = { "DoubleVortex" : DoubleVortex,     \
              "LinearXFlow2D" : LinearXFlow2D,     \
            }
    flowParams = [float(x) for x in lines[2].split("\t")[1:]]
    field = flows[flowName](*flowParams)

    # Get the initial and final positions
    p0 = np.array([float(x) for x in lines[3].split('\t')[1:]])
    pf = np.array([float(x) for x in lines[4].split('\t')[1:]])

    # Get the problem and solution parameters
    rhoE = float(lines[5].split("\t")[1])
    rhoT = float(lines[6].split("\t")[1])
    delJ = float(lines[7].split("\t")[1])
    initSpeedUbound = float(lines[8].split("\t")[1])
    compDist = float(lines[9].split("\t")[1])
    alpha = float(lines[10].split("\t")[1])

    # Create the solver and solve
    solver = AlphaFrontPropTESolver(field, p0, pf, rhoE, rhoT)

    solver.delJ = delJ
    solver.initSMax = initSpeedUbound
    solver.completionDist = compDist
    solver.alphaRad = alpha

    # The should be reasonable in pretty much all cases, so hardcode
    # them instead of setting for now.
    solver.nThInit = 50
    solver.nSInit = 50

    ## Initialize output file
    fileName = opts.o + "search_hist.txt"
    print("Initializing output file: " + fileName)
    outF = open(fileName, 'w')
    outF.write("Run name:\t%s Flow: %s\n"%(opts.o,field.name))
    outF.write("Run time:\t%s\n"%(time.strftime("%c")))
    outF.write("Run type:\t%s\n"%("2D-TIME-ENERGY"))
    outF.write("deltaJ:\t%s\n"%(solver.delJ))
    outF.write("Init P:\t%s\t%s\n"%(p0[0],p0[1]))
    outF.write("Final P:\t%s\t%s\n"%(pf[0],pf[1]))
    outF.write("\n")

    solver.GenInitPts()

    iters = 1000
    i = 0
    while(not solver.solutionFound and iters > 0):
        solver.RunJStep()
        iters = iters - 1

        # Write iteration to file solution
        print("Writing iteration " + str(i))
        outF.write("Iteration:\t%d\n"%(i))

        print("Writing points to file")
        points = solver.pointSets[-1]

        #np.set_printoptions(threshold='nan')
        #np.set_printoptions(linewidth=132)
        #print(points)

        outF.write("Num Pts:\t%d\n"%(len(points)))
        for iPt in range(len(points)):
            outF.write("%d\t%s\t%s\t%s\n"%(iPt,points[iPt,0],points[iPt,1],points[iPt,2]))

        print("Writing tris to file")
        tris = solver.surface
        outF.write("Num Tris:\t%d\n"%(len(tris)))
        for iTri in range(len(tris)):
            outF.write("%d\t%d\t%d\n"%(tris[iTri,0],tris[iTri,1],tris[iTri,2]))
        outF.write("\n")
        outF.flush()

        i = i + 1

    outF.close()

    ## Write out the solution path
    pathFilename = opts.o + 'optimal_trajectory.txt'
    print("Writing solution to " + pathFilename)
    pathFile = open(pathFilename, 'w')

    sol = solver.solution

    # Write out the header information
    pathFile.write("Final Cost:\t%s\n"%(str(solver.J)))
    pathFile.write("p0:\t"+"\t".join(str(x) for x in p0) +"\n")
    pathFile.write("pf:\t"+"\t".join(str(x) for x in pf) +"\n")
    pathFile.write("Num points:\t%d\n"%(len(sol)))

    # Write out each point
    for p in sol:
        pathFile.write("\t".join(str(x) for x in p) + "\n")
    pathFile.close()


if __name__ == "__main__":
    main(sys.argv)

