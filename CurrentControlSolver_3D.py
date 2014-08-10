# Solve a time optimal control problem in 3 dimensions and plot the result

# Nicholas Sharp - nsharp3@vt.edu

import os,sys,time
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as plt3D

from AlphaHullProp_3D import AlphaFrontPropTESolver
from FluidFuncs import *


def main(argV):

    ## Parse options
    usage = 'Run the 3D time-optimal alpha shape solver'

    parser = OptionParser(usage=usage)
    parser.add_option("", "-o", help="Output file name", metavar="string", default='3D_Time_Energy_soln.txt')

    (opts, args) = parser.parse_args(argV)

    ## Initial conditions
        
    p0 = np.array([0,0,0])
    pf = np.array([-1,1,1])

    # Background field
    field = LinearXFlow3D(.3)
    #field = Gyre3D(2, 2, 2, 2, 4)

    # Create the solver and solve
    solver = AlphaFrontPropTESolver(field, p0, pf)

    solver.delT = 0.15
    solver.sMax = 0.5

    solver.completionDist = 0.05
    scoopD = .25

    print("Scoop diameter = " + str(scoopD))
    #solver.alpha = 1.0/4.0 * scoopD**2
    solver.alphaRad = .5 * scoopD

    solver.nThInit = 20
    solver.nPhiInit = 20

    ## Initialize output file
    fileName = opts.o
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

