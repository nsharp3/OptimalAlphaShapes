# Solve a time optimal control problem in 3 dimensions and plot the result

# Nicholas Sharp - nsharp3@vt.edu

import os,sys,time
from optparse import OptionParser
import shutil

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as plt3D

from AlphaHullProp_3D import AlphaFrontPropTESolver
from FluidFuncs import *


def main(argV):

    ## Parse options
    usage = 'Run the 3D time-optimal alpha shape solver'

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
    delT:\tdelT
    speed:\tspeed
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

    print(lines)

    # Save a copy of the parameter file
    shutil.copyfile(fileName, opts.o + "input_params.txt")

    # Verify that this is the correct run type
    runType = lines[0].split("\t")[1]
    if runType != "3D-TIME":
        print("Incorrect run type, must be 3D-TIME")
        exit()

    # Get the type of flow and flow parameters
    flowName = lines[1].split("\t")[1]
    flows = { "Gyre3D" : Gyre3D     \
            }
    flowParams = [float(x) for x in lines[2].split("\t")[1:]]
    field = flows[flowName](*flowParams)

    # Get the initial and final positions
    p0 = np.array([float(x) for x in lines[3].split('\t')[1:]])
    pf = np.array([float(x) for x in lines[4].split('\t')[1:]])

    # Get the problem and solution parameters
    delT = float(lines[5].split("\t")[1])
    speed = float(lines[6].split("\t")[1])
    compDist = float(lines[7].split("\t")[1])
    alpha = float(lines[8].split("\t")[1])

    # Create the solver and solve
    solver = AlphaFrontPropTESolver(field, p0, pf)

    solver.delT = delT
    solver.sMax = speed
    solver.completionDist = compDist
    solver.alphaRad = alpha

    # The should be reasonable in pretty much all cases, so hardcode
    # them instead of setting for now.
    solver.nThInit = 50
    solver.nPhiInit = 50

    ## Initialize output file
    fileName = opts.o + "search_hist.txt"
    print("Initializing output file: " + fileName)
    outF = open(fileName, 'w')
    outF.write("Run name:\t%s Flow: %s\n"%(opts.o,field.name))
    outF.write("Run time:\t%s\n"%(time.strftime("%c")))
    outF.write("Run type:\t%s\n"%("3D-TIME"))
    outF.write("deltaT:\t%s\n"%(solver.delT))
    outF.write("Init P:\t%s\t%s\t%s\n"%(p0[0],p0[1],p0[2]))
    outF.write("Final P:\t%s\t%s\t%s\n"%(pf[0],pf[1],pf[2]))
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

        i = i + 1

    outF.close()

    ## Write out the solution path
    pathFilename = opts.o + 'optimal_trajectory.txt'
    print("Writing solution to " + pathFilename)
    pathFile = open(pathFilename, 'w')

    sol = solver.solution

    # Write out the header information
    pathFile.write("Final Time:\t%s\n"%(str(solver.t)))
    pathFile.write("p0:\t"+"\t".join(str(x) for x in p0) +"\n")
    pathFile.write("pf:\t"+"\t".join(str(x) for x in pf) +"\n")
    pathFile.write("Num points:\t%d\n"%(len(sol)))

    # Write out each point
    for p in sol:
        pathFile.write("\t".join(str(x) for x in p) + "\n")
    pathFile.close()


    '''
    # Plot solution

    fig = plt.figure(figsize=(8,8))
    fig.canvas.set_window_title("Surface View")
    ax = plt3D.Axes3D(fig)
    xLabel = ax.set_xlabel('X')
    yLabel = ax.set_ylabel('Y')
    zLabel = ax.set_zlabel('Z')
    points = solver.pointSets[-1]
    #ax.plot_trisurf(points[:,0],points[:,1],points[:,2], triangles=solver.surface, alpha=0.2)
    ax.scatter(points[:,0],points[:,1],points[:,2])
    plt.show()
    plt.close(fig)


    print(solver.solution)


    print("Plotting solution")
    fig = plt.figure(figsize=(14,14))
    fig.canvas.set_window_title("Solution Path")
    ax = plt3D.Axes3D(fig)
    xLabel = ax.set_xlabel('X (space)')
    yLabel = ax.set_ylabel('Y (space)')
    zLabel = ax.set_zlabel('t (time)')

    ax.set_xlim([-spaceLim,spaceLim])
    ax.set_ylim([-spaceLim,spaceLim])
    ax.set_zlim([-spaceLim,spaceLim])

    ax.plot(solver.solution[:,0], solver.solution[:,1], solver.solution[:,2], color ='red', linewidth=6)
    infoStr = 'steps = ' + str(solver.itNum) + '\ntime = ' + str(solver.t)
    ax.text2D(.05, .90, infoStr, transform=ax.transAxes)
    ax.scatter(p0[0],p0[1],p0[2], c='black', marker = 'o', s = 40)
    ax.scatter(pf[0],pf[1],pf[2], c='black', marker = 'o', s = 40)

    # Save several copies to create a long frame at the end of the movie
    while(i < solver.itNum + 20):
        plt.savefig(("%s/%06d.png"%(fileName,i)))
        i = i + 1

    plt.show()
    plt.close(fig)

    # Do a 2D plot of the solved system

    xArr = np.array(solver.solution[:,0])
    yArr = np.array(solver.solution[:,1])
    tArr = np.array(np.linspace(0, solver.t, len(solver.solution[:,0])))
    thArr = np.array(solver.solution[:,3])
    phiArr = np.array(solver.solution[:,4])

    fig2D = plt.figure(figsize=(14, 14))
    fig2D.canvas.set_window_title("Optimal Path Controls")

    thAx = plt.subplot2grid((2,4),(0,0),colspan=4)
    phiAx = plt.subplot2grid((2,4),(1,0),colspan=4)

    thAx.plot(tArr, thArr, color='black', linewidth=5)
    phiAx.plot(tArr, phiArr, color='black', linewidth=5)

    thAx.set_title("In-Plane Heading Control Function")
    thAx.set_ylabel("theta(t)")

    phiAx.set_title("Vertical Heading Control Function")
    phiAx.set_ylabel("phi(t)")
    phiAx.set_xlabel("Time")

    thAx.tick_params(axis='both')
    phiAx.tick_params(axis='both')

    while(i < solver.itNum + 40):
        plt.savefig(("%s/%06d.png"%(fileName,i)))
        i = i + 1

    plt.show()
    plt.close(fig2D)
    '''

if __name__ == "__main__":
    main(sys.argv)

