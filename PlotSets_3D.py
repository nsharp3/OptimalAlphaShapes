# Read in a solution file for a 3D front propagation run and create a plot image for each timestep

# Nicholas Sharp - nsharp3@vt.edu

import os,sys
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as plt3D


from FluidFuncs import *

'''
Solution file format:
(Currently assumes target is a static point)
Note that all data is TAB delimited
Comment lines start with '#'
===========================================
Run name:\tNAME
Run time:\tTIME
Run type:\t3D-TIME|2D-TIME-ENERGY
#Num iterations:\tNUM
Delta T:\tdelT
Init P:\tX\tY\tZ|X\tY
Final P:\tX\tY\tZ|X\tY

Iteration:\tNUM
Num Pts:\tNUM
IND\tX\tY\tZ
Num Tris:\tNUM
IND1\tIND2\tIND3
===========================================
'''

def main(argV):

    usage = 'Read in a solution file for 3D front propagation and create plot images'

    parser = OptionParser(usage=usage)
    parser.add_option("", "-o", help="prepend to outputs", metavar="string", default='')
    parser.add_option("", "-n", help="start file numbering at", type="int", default=0)
    parser.add_option("", "-s", help="skip to iteration", type="int", default=None)
    parser.add_option("", "--show-plots", help="Show the plots", action='store_true', default=False)
    parser.add_option("", "--draw-info", help="start numbering at", action='store_true', default=False)
    parser.add_option("", "--pdf", help="save an additional copy of the image as a pdf", action='store_true',default=False)
    parser.add_option("", "--make-movie", help="Run ffmpeg to make the move at the end", action='store_true', default=False)

    (opts, args) = parser.parse_args(argV)

    if len(args) != 2:
        print("Solution file must be given")
        exit()

    solFilename = args[1]
    print("Reading solution file from %s"%(solFilename))

    # Read through the file line by line
    # TODO: This means storing the entire file in memory, which could
    # get to be an issue for large runs
    solFile = open(solFilename,'r')
    # Get all lines without comments
    lines = []
    for line in solFile.readlines():
        if line[0] is not '#':
            lines.append(line)

    ## Read the header info
    runName = lines[0].strip().split('\t')[1]
    runTime = lines[1].strip().split('\t')[1]
    runType = lines[2].strip().split('\t')[1]
    #numIters = int(lines[3].strip().split('\t')[1])
    deltaT = float(lines[3].strip().split('\t')[1])

    if runType not in ['3D-TIME','2D-TIME-ENERGY']:
        print("Run type '%s' invalid"%(runType))
        exit()

    headerInitP = lines[4].strip().split('\t')[1:]
    initP = [float(x) for x in headerInitP]
    headerFinalP = lines[5].strip().split('\t')[1:]
    finalP = [float(x) for x in headerFinalP]

    # Print out the info from the header
    print("Header information:")
    print("\tRun name: %s"%(runName))
    print("\tRun time: %s"%(runTime))
    print("\tRun type: %s"%(runType))
    #print("\tNum iterations: %d"%(numIters))
    print("\tdeltaT: %s"%(str(deltaT)))
    print("\tInitial P: " + str(initP))
    print("\tFinal P: " + str(finalP))
    print("\n")

    ## Process the remaining iterations one at a time
    fileInd = 7
    iterInd = 0
    #for iterInd in range(numIters):
    while(True):

        print("\nPlotting iteration %d from file line %d"%(iterInd,fileInd))

        # First line contains the iteration number
        fileInd += 1
        if(fileInd >= len(lines)):
            print("Index off end of file, terminating")
            break

        ## Read in the points
        nPts = int(lines[fileInd].strip().split('\t')[1])
        print("\tReading %d points"%(nPts))
        fileInd += 1
        pts = []
        for iPt in range(nPts):
            lineData = [float(x) for x in lines[fileInd].strip().split('\t')]
            if lineData[0] != iPt:
                print("Points mis-numbered!")
                exit()
            pts.append(tuple(lineData[1:]))
            fileInd += 1


        ## Read in the triangles
        nTris = int(lines[fileInd].strip().split('\t')[1])
        print("\tReading %d triangles"%(nTris))
        fileInd += 1
        tris = []
        for iTri in range(nTris):
            lineData = [int(x) for x in lines[fileInd].strip().split('\t')]
            tris.append(tuple(lineData))
            fileInd += 1

        # Convert the points and triangles to numpy arrays
        pts = np.array(pts)
        tris = np.array(tris)

        # Possible skip plotting to satisfy the -s flag
        if iterInd < opts.s:
            print("\tSkipping plot iteration for -s flag")
            print("\tWill begin plotting at iteration %d"%(opts.s))
            iterInd += 1
            # Skip the blank line after the iteratoin
            fileInd += 1
            continue

        ## Make the plot

        # Plot solution
        # TODO: Get all of the text and line sizes right

        print("\tPlotting iteration")
        fig = plt.figure(figsize=(14,14))
        fig.canvas.set_window_title("Surface View")
        ax = plt3D.Axes3D(fig)
        xLabel = ax.set_xlabel('X')
        yLabel = ax.set_ylabel('Y')
        zLabel = ax.set_zlabel('Z')

        # TODO: Make the arguments to the program
        spaceLim = 3
        ax.set_xlim([-spaceLim,spaceLim])
        ax.set_ylim([-spaceLim,spaceLim])
        ax.set_zlim([-spaceLim,spaceLim])

        ax.plot_trisurf(pts[:,0],pts[:,1],pts[:,2], triangles=tris, color='red', shade=False)
        infoStr = 'iter = ' + str(iterInd) + '\nt = ' + str(iterInd*deltaT) + '\nnPts = ' + str(nPts) + '\nnTris = ' + str(nTris)
        if opts.draw_info:
            ax.text2D(.05, .90, infoStr, transform=ax.transAxes)


        # Save the plot as a png and optionall a pdf
        fileName = opts.o + "%06d"%(opts.n + iterInd)
        plt.savefig(fileName + ".png")
        if opts.pdf:
            plt.savefig(fileName + ".pdf")

        # Optionally show the plots
        if opts.show_plots:
            plt.show()

        plt.close(fig)

        iterInd += 1

        # Skip the blank line after the iteratoin
        fileInd += 1

    # Optionally render a movie
    if opts.make_movie:
        print("TODO: Implement this...")




if __name__ == "__main__":
    main(sys.argv)






