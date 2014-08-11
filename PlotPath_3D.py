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


def main(argV):

    # Set up the arguments to the program
    usage = 'Read in a solution file for 3D front propagation and create a 3D plot of the optimal trajectory as well as 2D plots of the optimal controls'
    parser = OptionParser(usage=usage)

    parser.add_option("", "--format", help="output file format", metavar="string", default='pdf')
    parser.add_option("", "--space-lim", help="spatial limits for the plot", type="float", default=3.0)
    parser.add_option("", "--fig-size", help="set the size of the plotting window", type="int", default=6)
    parser.add_option("", "--make-movie", help="Run ffmpeg to make a movie at the end (with default settings that are usually decent)", action='store_true', default=False)

    solDir = args[1]

    # Read in the solution data
    solFilename = soldDir + 'search_hist.txt'
    print("Reading solution file from %s"%(solFilename))
    # Get all lines without comments
    lines = []
    for line in solFile.readlines():
        if line[0] is not '#':
            lines.append(line.strip())

    # Read the initial and final points
    p0 = [float(x) for x in lines[0].split('\t')[1:]]
    pf = [float(x) for x in lines[1].split('\t')[1:]]

    # Read the actual points in the trajectory
    nPts = int(lines[2].split('\t')[1])
    sol = []
    for i in range(3,3+nPts):
        sol.append([float(x) for x in lines[i].split('\t')])
    sol = np.array(sol)

    # Make a 3D plot of the solution
    print("Making 3D plot of path")
    fig = plt.figure(figsize=(14,14))
    fig.canvas.set_window_title("Solution Path")
    ax = plt3D.Axes3D(fig)
    xLabel = ax.set_xlabel('X (space)')
    yLabel = ax.set_ylabel('Y (space)')
    zLabel = ax.set_zlabel('t (time)')

    spaceLim = opts.space_lim
    ax.set_xlim([-spaceLim,spaceLim])
    ax.set_ylim([-spaceLim,spaceLim])
    ax.set_zlim([-spaceLim,spaceLim])

    # TODO: Add an info string with num iterations, final time, etc
    ax.plot(sol[:,0], sol[:,1], sol[:,2], color ='red', linewidth=6)
    ax.scatter(p0[0],p0[1],p0[2], c='black', marker = 'o', s = 40)
    ax.scatter(pf[0],pf[1],pf[2], c='black', marker = 'o', s = 40)

    # Save the plot
    fileName = solDir + 'pointcloud' + "%06d"%(opts.n + iterInd) + '.' + opts.format
    plt.savefig(fileName)


    if opts.make_movie:
        print("TODO: Implement make-movie")


if __name__ == "__main__":
    main(sys.argv)
