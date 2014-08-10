# Solve time-energy current optimal control problems using alpha-shapes
# to manage fronts.

# Nicholas Sharp - nsharp@vt.edu

import numpy as np

import getopt, sys, os, stat, subprocess, string

from FluidFuncs import *


class AlphaFrontPropTESolver:

    # Initialize the solver
    def __init__(self, fluidFunc, p0, pf):
        
        self.p0 = p0
        self.pf = pf
        self.fluidFunc = fluidFunc
       
        self.t = 0

        # Initialization parameters
        self.nThInit = 50
        self.nPhiInit = 50
        self.sMax = 3
        
        
        # Algorithm parameters
        self.solutionFound = False
        self.itNum = 0
        self.delT = 0.1
        self.delD = .05

        # CGAL documentation appears inconsistent with Edelsbrunner papers. Edelsbrunner uses 
        # alpha = radius of scoop (see formal mathematical definition). CGAL uses alpha = squared
        # radius of scoop. This value is the Edelsbrunner value, conversion below.
        self.alphaRad = -1

        # The distance from ayn
        self.completionDist = -1

        # Stored as a python list (for efficent appends) of numpy Nx6 arrays 
        # representing each iteration. Each of the N points is stored as
        # [x ,y, t, th, s, prevIndex] where the first 5 elements are reals
        # and the last is an index into the array from the previous iteration
        self.pointSets = []
        
        # The suface corresponding to the current point set, computed by the 
        # alpha-shape algorithm. Stored an Nx3 array, where N is the number
        # of triangles and each triangle is stored as [pt0, pt1, pt2], which are
        # indices in to the current point set.
        self.surface = np.array([])
    
        # Internal members
        self.interpStart = -1

    # Generate the points to intialize the search
    def GenInitPts(self):
    
        nInitPts = self.nThInit * self.nPhiInit
        initPts = np.empty([nInitPts + 1,6])
        
        thArr = np.linspace(0, 2*np.pi, self.nThInit, endpoint=False)
        phiArr = np.linspace(0, np.pi, self.nPhiInit+1, endpoint=False)
        
        # Note the indexing tricks to exclude the start and ending point
        # of phiArr

        index = 0
        for iTh in range(self.nThInit):
            for iS in range(1,self.nPhiInit+1):
                initPts[index] = np.array([self.p0[0], self.p0[1], self.p0[2], thArr[iTh], phiArr[iS], -1])
                index = index + 1
                
        self.pointSets.append(initPts)
    
    def RunJStep(self):
    
        self.itNum = self.itNum + 1
        print("\n===== Beginning search iteration " + str(self.itNum) + " =====")
        
        # Interpolate points along the surface
        print("Interpolating points along surface")
        self.InterpolateSurface()
    
        # Propagate points
        print("Propagating front points")
        self.PropagatePoints()
        
        # Compute an alpha-hull 
        print("Computing alpha hull")
        self.ComputeHull()
        
        self.t = self.t + self.delT

        # Check if the set containst the target
        print("Checking for completion")
        self.CheckCompletion()
       
        # Generate the solution set, if needed 
        if(self.solutionFound):
            print("SOLUTION FOUND!!!")
            self. ReconstructSolution()
    
    def PropagatePoints(self):
        
        currPts = self.pointSets[-1]
      
        newPts = []

        # Count how many times each point appears in at least one triangle of the alpha-shape
        counts = np.zeros(len(currPts))
        counts[self.interpStart:] = 1 # The newly interpolated points
        for tri in self.surface:
            counts[tri[0]] = counts[tri[0]] + 1
            counts[tri[1]] = counts[tri[1]] + 1
            counts[tri[2]] = counts[tri[2]] + 1

        for i in range(len(currPts)):
            if(counts[i] > 0 or self.itNum == 1):
                deltas = self.OptimalDeltas(currPts[i])
                newPt = currPts[i] + deltas
                newPt[5] = i
                newPts.append(newPt)
            
        print("   " + str(len(newPts)) + " points were propagated, " + str(len(currPts) - len(newPts)) + " were not")
        self.pointSets.append(np.array(newPts))
        
    # Derived via Optimal Control Theory. See associated writeup.

    def OptimalDeltas(self, X):
    
        th = X[3]
        phi = X[4]

        xDot = self.fluidFunc.Ux(X[0], X[1], X[2], self.t) + self.sMax*np.cos(th)*np.sin(phi)
        yDot = self.fluidFunc.Uy(X[0], X[1], X[2], self.t) + self.sMax*np.sin(th)*np.sin(phi)
        zDot = self.fluidFunc.Uz(X[0], X[1], X[2], self.t) + self.sMax*np.cos(phi)

        # Get the derivatives all at once
        dUxdx = self.fluidFunc.dUxdx(X[0],X[1],X[2],self.t)
        dUxdy = self.fluidFunc.dUxdy(X[0],X[1],X[2],self.t)
        dUxdz = self.fluidFunc.dUxdz(X[0],X[1],X[2],self.t)
        
        dUydx = self.fluidFunc.dUydx(X[0],X[1],X[2],self.t)
        dUydy = self.fluidFunc.dUydy(X[0],X[1],X[2],self.t)
        dUydz = self.fluidFunc.dUydz(X[0],X[1],X[2],self.t)
        
        dUzdx = self.fluidFunc.dUzdx(X[0],X[1],X[2],self.t)
        dUzdy = self.fluidFunc.dUzdy(X[0],X[1],X[2],self.t)
        dUzdz = self.fluidFunc.dUzdz(X[0],X[1],X[2],self.t)

        # See VSGC 2014 paper

        if (phi < 0.001):
            # Avoid numerical blowup for degenerate phi
            thDot = 0
        else:
            thDot = 0.5*(   (dUxdx - dUydy)*np.sin(2*th)
                            - (dUxdy + dUydx)*np.cos(2*th) 
                            - dUxdy + dUydx
                            + 2*(1/np.tan(phi))*(dUzdx*np.sin(th) - dUzdy*np.cos(th))
                        )

        phiDot = 0.25*( - np.sin(2*phi)*( (dUxdx - dUydy)*np.cos(2*th)
                                         + dUxdx + dUydy - 2*dUzdz
                                         - (dUxdy + dUydx)*np.sin(2*th)
                                         )
                        -2*np.cos(th)* ( (dUxdz + dUzdx)*np.cos(2*phi)
                                         - dUxdz + dUzdx
                                       )
                        -2*np.sin(th)* ( (dUydz + dUzdy)*np.cos(2*phi)
                                         - dUydz + dUzdy
                                       )
                        ) 


        delX = xDot * self.delT
        delY = yDot * self.delT
        delZ = zDot * self.delT
        delTh = thDot * self.delT
        delPhi = phiDot * self.delT

        return np.array([delX, delY, delZ, delTh, delPhi, 0])
   
    # Wraps CGAL's Alpha-Hull implmementation
    def ComputeHull(self):
        
        hullProc = subprocess.Popen(os.path.join(os.path.dirname(__file__),"CGAL_Alpha_Wrapper"), 
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        points = self.pointSets[-1]
        
        # See above for inconsistency of definitions
        CGAL_Alpha = self.alphaRad**2

        # Efficiently (ish) build the giant string to pass to stdin
        inStr = str(CGAL_Alpha) + "\n" + str(len(points)) + "\n"  + ''.join([str(point[0]) + " " + str(point[1]) + " "  + str(point[2]) + "\n" for point in points])
        
        # Let the algorithm run. There will be extraneous output
        (outStr, errStr) = hullProc.communicate(input=inStr)
        hullProc.stdin.close()
   
        if(len(errStr) > 0):
            print("   Error output from CGAL Alpha Hull:")
            print("      " + errStr)

        # Parse output to surface
        triList = []
        for line in string.split(outStr, '\n'):
            splitStr = line.split(" ")
            if (len(splitStr) == 3):
                nums = map(int, line.split(" "))
                triList.append([nums[0],nums[1],nums[2]])

        self.surface = np.array(triList)
       
        # TODO: Walk graph to remove interior faces

    def InterpolateSurface(self):
    
        # For every edge (in every triangle), if it is longer than some threshold place
        # a new point at its center

        points = self.pointSets[-1]
        
        thresh = self.alphaRad      # half of the diameter of the ball
        print("   Interpolation threshold = " + str(thresh))
        
        newPoints = []

        # Since the test happens on an edge-by-edge basis, each edge could get interpolated twice. Use 
        # a set to prevent this
        interpSet = set()

        for tri in self.surface:
            for i in (0,1,2):
                p1 = points[tri[i]]
                p2 = points[tri[(i+1)%3]]
    
                if((p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]) in interpSet or (p2[0],p2[1],p2[2],p1[0],p1[1],p1[2]) in interpSet):
                    continue

                sLen = np.linalg.norm((p1-p2)[0:3])
                if (sLen > thresh):
                    newP = (p1 + p2) / 2
                    # FIXME: Need to calculate phi properly 
                    # Mean of circular quantities:
                    meanTh = np.arctan2(.5*(np.sin(p1[3])+np.sin(p2[3])),.5*(np.cos(p1[3])+np.cos(p2[3])))
                    prev = p1[5] # p1, arbitrary
                    newP[3] = meanTh
                    newP[5] = prev
                    newPoints.append(newP)
       
                interpSet.add((p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]))

        print("Interpolation added " + str(len(newPoints)) + " new points.")
       
        self.interpStart = len(points)

        if ( len(newPoints) > 0):
            self.pointSets[-1] = np.vstack([points,np.array(newPoints)])

    
    def LineVecCross(self, pFrom, pTo, pTest):
        vTest0 = pTest[0] - pFrom[0]
        vTest1 = pTest[1] - pFrom[1]
        vTo0 = pTo[0] - pFrom[0]
        vTo1 = pTo[1] - pFrom[1]

        return vTest0 * vTo1 - vTest1 * vTo0

    # Checks if the final point is contained in the x-y projection of the reachable set
    def CheckCompletion(self):

        # Check if the point is sufficiently close to any vertex. Assumption: this means
        # the reachable set is about to include that point

        # TODO: This is a bad way of doing things. Find something better        
        # Find the closest point in the current pointset
        
        minDist = (np.apply_along_axis(np.linalg.norm,1,(self.pointSets[-1][:,0:3]-self.pf))).min()
        
        print("    Min dist = "+str(minDist))    

        self.solutionFound =  minDist <= self.completionDist


    def ReconstructSolution(self):

        if(not self.solutionFound):
            print("Reconstruct solution called wrong!!")

        pointTrail = []

        # Find the closest point in the current pointset
        minInd = (np.apply_along_axis(np.linalg.norm,1,(self.pointSets[-1][:,0:3]-self.pf))).argmin()
       
        print("Closest point:")
        print(self.pointSets[-1][minInd])

        pointTrail.append(self.pointSets[-1][minInd])
        nextInd = self.pointSets[-1][minInd,5]

        # Walk the list backwards
        i = -2
        while(nextInd != -1):
            pointTrail.append(self.pointSets[i][nextInd])
            nextInd = self.pointSets[i][nextInd,5]
            i = i - 1
        
        pointTrail.reverse()

        self.solution = np.array(pointTrail)

