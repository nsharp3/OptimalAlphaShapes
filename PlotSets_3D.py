# Read in a search history file for a 3D front propagation run and create a plot image for each timestep

# Nicholas Sharp - nsharp3@vt.edu

import os,sys
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as plt3D


from FluidFuncs import *

'''
Search history file format:
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
Speed|Max Speed:\tspeed

Iteration:\tNUM
Num Pts:\tNUM
IND\tX\tY\tZ
Num Tris:\tNUM
IND1\tIND2\tIND3
===========================================
'''

# RC params to make pretty plots
# See: http://damon-is-a-geek.com/publication-ready-the-first-time-beautiful-reproducible-plots-with-matplotlib.html 
from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
from matplotlib.ticker import MaxNLocator
my_locator = MaxNLocator(6)

def main(argV):

    # Set up the arguments to the program
    usage = '''python PlotSets.py SOLUTION_DIRECTORY\n
Read in a solution file for 3D front propagation and create plot images. By default, every iterations is plotted and saved. Use -s to plot a singler iteration.

The SOLUTION_DIRECTOY should be the directory containing the data files for the run. Plots will also be saved there.'''

    parser = OptionParser(usage=usage)
    parser.add_option("", "-n", help="start file numbering at", type="int", default=0)
    parser.add_option("", "-s", help="skip to iteration", type="int", default=None)
    parser.add_option("", "--format", help="output file format", metavar="string", default='png')
    parser.add_option("", "--space-lim", help="spatial limits for the plot", type="float", default=3.0)
    parser.add_option("", "--fig-size", help="set the size of the plotting window", type="int", default=6)
    parser.add_option("", "--show-plots", help="Show the plots", action='store_true', default=False)
    parser.add_option("", "--draw-info", help="start numbering at", action='store_true', default=False)
    parser.add_option("", "--make-movie", help="Run ffmpeg to make a movie at the end (with default settings that are usually decent)", action='store_true', default=False)
    parser.add_option("", "--rep-last-frame", help="Produce N copies of the last frame to give a pausing effect.", type="int", default=1)
    parser.add_option("", "--point-cloud", help="Also make plots of the unconnected point cloud", action='store_true', default=False)

    (opts, args) = parser.parse_args(argV)

    # Validate input
    if len(args) != 2:
        print("Solution file must be given")
        exit()
    if opts.format not in ['png','pdf','eps']:
        print("format unknown")
        exit()

    solDir = args[1]

    solFilename = solDir + 'search_hist.txt'
    print("Reading solution file from %s"%(solFilename))

    # Program parameters
    surfColor = '#3366FF'

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
    delta = float(lines[3].strip().split('\t')[1])

    if runType not in ['3D-TIME','2D-TIME-ENERGY']:
        print("Run type '%s' invalid"%(runType))
        exit()

    # Don't do anything with the speed for now

    headerInitP = lines[5].strip().split('\t')[1:]
    initP = [float(x) for x in headerInitP]
    headerFinalP = lines[6].strip().split('\t')[1:]
    finalP = [float(x) for x in headerFinalP]

    # Print out the info from the header
    print("Header information:")
    print("\tRun name: %s"%(runName))
    print("\tRun time: %s"%(runTime))
    print("\tRun type: %s"%(runType))
    #print("\tNum iterations: %d"%(numIters))
    print("\tdelta: %s"%(str(delta)))
    print("\tInitial P: " + str(initP))
    print("\tFinal P: " + str(finalP))
    print("\n")

    ## Process the remaining iterations one at a time
    fileInd = 8
    iterInd = 0
    #for iterInd in range(numIters):
    while(True):
        

        print("\nPlotting iteration %d from file line %d"%(iterInd,fileInd))

        # First line contains the iteration number
        fileInd += 1
        if(fileInd >= len(lines)):
            print("Index off end of file, terminating")
            break

        if opts.rep_last_frame > 1 and iterInd > 0:
            plt.close(fig)

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
        if opts.s != None and iterInd > opts.s:
            print("\tIteration is after the one requested from -s, exiting")
            break

        ## Make the plot

        # Plot solution
        # TODO: Get all of the text and line sizes right

        print("\tPlotting iteration")
        figSize = (opts.fig_size,opts.fig_size)
        fig = plt.figure(figsize=figSize)
        fig.canvas.set_window_title("Surface View")
        ax = plt3D.Axes3D(fig)
        xLabel = ax.set_xlabel('X')
        yLabel = ax.set_ylabel('Y')
        if runType == "3D-TIME":
            zLabel = ax.set_zlabel('Z')
        if runType == "2D-TIME-ENERGY":
            zLabel = ax.set_zlabel('T')
        #ax.xaxis.set_major_locator(my_locator)
        #ax.yaxis.set_major_locator(my_locator)
        #ax.zaxis.set_major_locator(my_locator)

        # TODO: Make these arguments to the program or something
        # (also used below in point_cloud)
        spaceLim = opts.space_lim
        ax.set_xlim([-spaceLim,spaceLim])
        ax.set_ylim([-spaceLim,spaceLim])
        if runType == "3D-TIME":
            ax.set_zlim([-spaceLim,spaceLim])
        if runType == "2D-TIME-ENERGY":
            ax.set_zlim([0,spaceLim])

        # Draw the actual plot
        ax.plot_trisurf(pts[:,0],pts[:,1],pts[:,2], triangles=tris, color=surfColor, shade=False, alpha=0.5, linewidth=0.1)

        # Draw the target
        if runType == "2D-TIME-ENERGY":
            ax.plot([finalP[0],finalP[0]],[finalP[1],finalP[1]],[0,spaceLim], c='black')
        else:
            ax.scatter([finalP[0]],[finalP[1]], [finalP[2]], color='black', s=70, marker='*')

        if runType == "2D-TIME-ENERGY":
            costWord = "J"
        if runType == "3D-TIME":
            costWord = "t"

        infoStr = 'iter = ' + str(iterInd) + '\n '+ costWord +' = ' + str(iterInd*delta) + '\nnPts = ' + str(nPts) + '\nnTris = ' + str(nTris)
        if opts.draw_info:
            ax.text2D(.05, .85, infoStr, transform=ax.transAxes, bbox=dict(facecolor='grey', alpha=0.5))


        # Save the plot as the appropriate file type
        fileName = solDir + "%06d"%(opts.n + iterInd) + '.' + opts.format
        print("\tSaving plot file " + fileName)
        plt.savefig(fileName)

        # Optionally show the plots
        if opts.show_plots:
            plt.show()

        # If this is the last frame, we made need to leave it around so
        # we can create duplicate plots of it at the end.
        if not opts.rep_last_frame > 0:
            plt.close(fig)

        # Optionally make an unconnected point cloud plot
        if opts.point_cloud:

            print("\tPlotting point cloud")
            fig = plt.figure(figsize=figSize)
            fig.canvas.set_window_title("Point Cloud")
            ax = plt3D.Axes3D(fig)
            xLabel = ax.set_xlabel('X')
            yLabel = ax.set_ylabel('Y')
            zLabel = ax.set_zlabel('Z')

            # TODO: Make these arguments to the program or something
            ax.set_xlim([-spaceLim,spaceLim])
            ax.set_ylim([-spaceLim,spaceLim])
            ax.set_zlim([-spaceLim,spaceLim])
            ax.xaxis.set_major_locator(my_locator)
            ax.yaxis.set_major_locator(my_locator)
            ax.zaxis.set_major_locator(my_locator)

            # Draw the actual plot
            ax.scatter(pts[:,0],pts[:,1],pts[:,2], color='black', s=4)
            infoStr = 'iter = ' + str(iterInd) + '\nt = ' + str(iterInd*delta) + '\nnPts = ' + str(nPts) + '\nnTris = ' + str(nTris)
            if opts.draw_info:
                ax.text2D(.05, .85, infoStr, transform=ax.transAxes, bbox=dict(facecolor='grey', alpha=0.5))

            # Save the plot
            fileName = solDir + 'pointcloud' + "%06d"%(opts.n + iterInd) + '.' + opts.format
            plt.savefig(fileName)

            # Optionally show the plots
            if opts.show_plots:
                plt.show()

            plt.close(fig)


        # Skip the blank line after the iteratoin
        fileInd += 1

        iterInd += 1

    # Create duplicates of the last frame
    if opts.rep_last_frame > 1:
        for i in range(opts.rep_last_frame):
            print("\tPlotting last frame")
            fileName = solDir + "%06d"%(opts.n + iterInd + i) + '.' + opts.format
            print("\tSaving plot file " + fileName)
            plt.savefig(fileName)


    # Optionally render a movie
    if opts.make_movie:
        print("TODO: Implement this...")

if __name__ == "__main__":
    main(sys.argv)






