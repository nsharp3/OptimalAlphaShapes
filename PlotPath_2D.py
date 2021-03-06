# Read in a solution file for a 3D front propagation run and plot the
# optimal trajectory

# Nicholas Sharp - nsharp3@vt.edu

import os,sys
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as plt3D


from FluidFuncs import *

'''
Solution file format
p0:\tx\ty\tz
pf:\tx\ty\tz
Num points:\tNUM
Xi\tYi\tZi\tTHi\tPHIi
...
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

# Pleasant colorscheme:
niceBlue = '#3366FF'
niceRed  = '#CC0000'

def main(argV):

    # Set up the arguments to the program
    usage = '''python PlotSets.py SOLUTION_DIRECTORY\n
Read in a solution file for 2D time-energy front propagation and create a 3D plot of the optimal trajectory as well as 2D plots of the optimal controls

The SOLUTION_DIRECTOY should be the directory containing the data files for the run. Plots will also be saved there.'''
    parser = OptionParser(usage=usage)

    parser.add_option("", "--format", help="output file format", metavar="string", default='pdf')
    parser.add_option("", "--space-lim", help="spatial limits for the plot", type="float", default=3.0)
    parser.add_option("", "--fig-size", help="set the size of the plotting window", type="int", default=6)
    parser.add_option("", "--show-plots", help="Show the plots", action='store_true', default=False)
    parser.add_option("", "--make-movie", help="Run ffmpeg to make a movie at the end (with default settings that are usually decent)", action='store_true', default=False)

    (opts, args) = parser.parse_args(argV)

    solDir = args[1]

    # Read in the solution data
    solFilename = solDir + 'optimal_trajectory.txt'
    print("Reading solution file from %s"%(solFilename))
    solFile = open(solFilename,'r')
    # Get all lines without comments
    lines = []
    for line in solFile.readlines():
        if line[0] is not '#':
            lines.append(line.strip())

    # Read the initial and final points
    finalT = float(lines[0].split('\t')[1])
    p0 = [float(x) for x in lines[1].split('\t')[1:]]
    pf = [float(x) for x in lines[2].split('\t')[1:]]

    # Read the actual points in the trajectory
    nPts = int(lines[3].split('\t')[1])
    sol = []
    for i in range(4,4+nPts):
        sol.append([float(x) for x in lines[i].split('\t')])
    sol = np.array(sol)

    # Make a 3D plot of the solution
    print("Making 3D plot of path")
    figSize = (6,6)
    fig = plt.figure(figsize=figSize)
    fig.canvas.set_window_title("Solution Path")
    ax = plt3D.Axes3D(fig)
    xLabel = ax.set_xlabel('X')
    yLabel = ax.set_ylabel('Y')
    zLabel = ax.set_zlabel('Z')

    spaceLim = opts.space_lim
    ax.set_xlim([-spaceLim,spaceLim])
    ax.set_ylim([-spaceLim,spaceLim])
    ax.set_zlim([0,spaceLim])

    # TODO: Add an info string with num iterations, final time, etc
    # The solution path
    ax.plot(sol[:,0], sol[:,1], sol[:,2], color=niceRed, linewidth=6)
    # The target line
    ax.plot([pf[0],pf[0]],[pf[1],pf[1]],[0,spaceLim], c='black')

    # Save the plot
    fileName = solDir + 'optimal_traj' + '.' + opts.format
    plt.savefig(fileName)

    if opts.show_plots:
        plt.show()
    plt.close(fig)

    if opts.make_movie:
        print("TODO: Implement make-movie")

    # Assemble the data
    xArr = np.array(sol[:,0])
    yArr = np.array(sol[:,1])
    tArr = np.array(sol[:,2])
    thArr = np.array(sol[:,3])
    sArr = np.array(sol[:,4])

    ## Make 2D plot of the path with the current in the background
    print("TODO")


    ## Make 2D plots of the optimal controls

    # Create the figure and axes
    fig2D = plt.figure(figsize=figSize)
    fig2D.canvas.set_window_title("Optimal Path Controls")
    thAx = plt.subplot2grid((2,4),(0,0),colspan=4)
    sAx = plt.subplot2grid((2,4),(1,0),colspan=4)

    # Plot the lines
    thAx.plot(tArr, thArr, color=niceBlue, linewidth=5)
    sAx.plot(tArr, sArr, color=niceBlue, linewidth=5)

    # Titles and labels
    thAx.set_title("Heading Control Function")
    thAx.set_ylabel("$\\theta(t)$")

    sAx.set_title("Speed Control Function")
    sAx.set_ylabel("$s(t)$")
    sAx.set_xlabel("Time")

    thAx.tick_params(axis='both')
    sAx.tick_params(axis='both')

    # Save the plot
    fileName = solDir + 'optimal_controls' + '.' + opts.format
    plt.savefig(fileName)

    if opts.show_plots:
        plt.show()
    plt.close(fig2D)



if __name__ == "__main__":
    main(sys.argv)






