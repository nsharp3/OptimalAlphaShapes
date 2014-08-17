# Solve time-energy current optimal control problems using alpha-shapes
# to manage fronts.

# Nicholas Sharp - nsharp@vt.edu

import numpy as np
np.set_printoptions(threshold=np.nan)
np.set_printoptions(linewidth=200)

import getopt, sys, os, stat, subprocess, string
from math import sin,cos,tan
import random as rand

from FluidFuncs import *


class AlphaFrontPropTESolver:

    # Initialize the solver
    def __init__(self, fluidFunc, p0, pf, rhoE, rhoT):
        
        self.p0 = p0
        self.pf = pf
        self.rhoE = rhoE
        self.rhoT = rhoT
        self.fluidFunc = fluidFunc
       
        self.J = 0

        # Initialization parameters
        self.nThInit = 50
        self.nSInit = 50
        self.initSMax = 5 
        
        # Algorithm parameters
        self.solutionFound = False
        self.itNum = 0
        self.delJ = 0.1
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
    
        initPts = []
        
        thArr = np.linspace(0, 2*np.pi, self.nThInit, endpoint=False)
        sArr = np.linspace(0, self.initSMax, self.nSInit, endpoint=False)
        
        for iTh in range(self.nThInit):
            for iS in range(self.nSInit):
                initPts.append((self.p0[0], self.p0[1], 0, thArr[iTh], sArr[iS], -1))

        initPts = np.array(initPts)
                
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
        
        self.J = self.J + self.delJ

        # Compute an alpha shape 
        print("Computing alpha shape")
        self.ComputeHull()
        
        # Reconcile the shape to hull
        print("Reconciling alpha shape to hull")
        hullFound = self.ReconcileHull()
        if !hullFound:
            print("Warning: alpha shape could not be reconciled to a hull. Exiting.")
            exit()

        # Check if the set containst the target
        print("Checking for completion")
        self.CheckCompletion()
       
        # Generate the solution set, if needed 
        if(self.solutionFound):
            print("SOLUTION FOUND!!!")
            self.ReconstructSolution()
    
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

                # As mentioned in the paper, this is not the
                # time-optimal special case, so we must retain the
                # points from previous time steps if they still fall on
                # the reachable set.
                newPt = np.copy(currPts[i])
                newPt[5] = i
                newPts.append(newPt)
  
        np.set_printoptions(threshold=np.nan)
        np.set_printoptions(linewidth=200)
        #print(np.array(newPts))
        #exit()

        print("   " + str(len(newPts)) + " points were propagated, " + str(2*len(currPts) - len(newPts)) + " were not")
        self.pointSets.append(np.array(newPts))
        
    # Derived via Optimal Control Theory. See associated writeup.

    def OptimalDeltas(self, X):
    
        t = X[2]
        th = X[3]
        s = X[4]

        xDot = self.fluidFunc.Ux(X[0], X[1], t) + s*cos(th)
        yDot = self.fluidFunc.Uy(X[0], X[1], t) + s*sin(th)

        # Get the derivatives all at once
        dUxdx = self.fluidFunc.dUxdx(X[0], X[1], t)
        dUxdy = self.fluidFunc.dUxdy(X[0], X[1], t)
        
        dUydx = self.fluidFunc.dUydx(X[0], X[1], t)
        dUydy = self.fluidFunc.dUydy(X[0], X[1], t)

        # Compute optimal control changes
        # See paper explanation (or Mathematica notebook TE_Derivations.nb)
        thDot = -dUxdy*cos(th)**2 + dUydx*sin(th)**2            \
                    +   (dUxdx-dUydy)*cos(th)*sin(th)
        
        sDot = - 0.5 * s * (dUxdx*cos(th)**2 + dUydy*sin(th)**2            \
                                +   (dUxdy+dUydx)*cos(th)*sin(th))

        delT = self.delJ / (self.rhoE * s**3 + self.rhoT) 


        delX = xDot * delT
        delY = yDot * delT
        delTh = thDot * delT
        delS = sDot * delT

        return np.array([delX, delY, delT, delTh, delS, 0])
   
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
  
        #print("Output from CGAL:")
        #print(outStr)
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
       
    # Given an alpha shape, reconcile it to an alpha hull or report that
    # it cannot be done. See AlphaReconciliation.txt for a full
    # explanation.
    def ReconcileHull(self):

        pts = self.pointSets[-1]
        aShape = self.surface

        # Compute bounds
        mins = np.min(pts,0)
        maxs = np.max(pts,0)

        # Compute p0
        r1 = rand.random()
        r2 = rand.random()
        p0x = mins[0] - 0.5 * (maxs[0] - mins[0])
        p0y = mins[1] - r1 * (maxs[1] - mins[1])
        p0z = mins[2] - r2 * (maxs[2] - mins[2])

        # Compute p1
        rInd = rand.randint(0,len(aShape)-1)
        p1x = 0.0
        p1y = 0.0
        p1z = 0.0
        for i in range(3):
            p1x += pts[aShape[rInd,i],0]
            p1y += pts[aShape[rInd,i],1]
            p1z += pts[aShape[rInd,i],2]
        p1x = p1x / 3.0
        p1y = p1y / 3.0
        p1z = p1z / 3.0

        # Adjust p1 so that it definitely intersects the triangle it was
        # generated from
        p1x += 0.1 * (p1x - p0x)
        p1y += 0.1 * (p1y - p0y)
        p1z += 0.1 * (p1z - p0z)

        # Find the triangular facet which intersects the line nearest p0
        

    # Tests if a ray intersects a triangle. Returns none if there is no
    # intersection, or float s >= 0 which is a paramaterization of the
    # intersection point along p0 --> p1
    # 
    # p0 and p1 are points. tri should be a 3x3 array of point x coord
    #
    # From http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()
    # TODO do they have something I can cite?
    def RayIntersectTriangle(p0, p1, tri):
        
        EPS = 0.000000001

        # Get triangle edge vectors and plane normal
        u = tri[1,:] - tri[0,:]
        v = tri[2,:] - tri[0,:]
        n = np.cross(u,v)

        d = p1 - p0     # Ray direction vector
        w0 = p0 - tri[0,:]
        a = -np.dot(n,w0)
        b = np.dot(n,d)

        # The ray is inplane or disjoint from plane
        if np.abs(b) < EPS:
            return None

        r = a / b       # The parameterization of the intersection point along the ray
        
        # Test if the intersection point is "behind" the ray
        if r < 0:
            return None

        intPoint = p0 + r*d

        # Test if the intersection point is inside the triangle
        uu = np.dot(u,u)
        uv = np.dot(u,v)
        vv = np.dot(v,v)
        w = intPoint - tri[0,:]
        wu = np.dot(w,u)
        wv = np.dot(w,v)
        D = uv*uv - uu*vv

        # Compute and test the parametric coordinates
        s = (uv * wv - vv * wu) / D
        if s < 0 or s > 1:
            return None
        t = (uv * wu - uu * wv) / D
        if t < 0 or (s+t) > 1:
            return None

        # This intersection is valid! Return the parameterization
        return r


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

        # Project every triangle onto the x-y plane and check if the point
        # is within that triangle
        points = self.pointSets[-1]
        for tri in self.surface:

            # Cross products
            p1 = points[tri[0]]
            p2 = points[tri[1]]
            p3 = points[tri[2]]
            c0 = self.LineVecCross(p1, p2, self.pf) 
            c1 = self.LineVecCross(p2, p3, self.pf) 
            c2 = self.LineVecCross(p3, p1, self.pf) 

            # Check signs of cross products
            if ((c0 > 0 and c1 > 0 and c2 > 0) or (c0 < 0 and c1 < 0 and c2 < 0)):
                self.solutionFound = True;
                return

    def ReconstructSolution(self):

        if(not self.solutionFound):
            print("Reconstruct solution called wrong!!")

        pointTrail = []

        # Find the closest point in the current pointset
        minInd = (np.apply_along_axis(np.linalg.norm,1,(self.pointSets[-1][:,0:2]-self.pf))).argmin()
       
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

