# Solve time-energy current optimal control problems using alpha-shapes
# to manage fronts.

# Nicholas Sharp - nsharp@vt.edu

import numpy as np

import getopt, sys, os, stat, subprocess, string
from math import sin,cos,tan
import random as rand
import Queue
from collections import defaultdict

from FluidFuncs import *


class AlphaFrontPropTSolver:

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

        initPts = []

        thArr = np.linspace(0, 2*np.pi, self.nThInit, endpoint=False)
        phiArr = np.linspace(0, np.pi, self.nPhiInit+1, endpoint=False)

        # Note the indexing tricks to exclude the start and ending point
        # of phiArr

        for iTh in range(self.nThInit):
            for iS in range(1,self.nPhiInit+1):
                initPts.append((self.p0[0], self.p0[1], self.p0[2], thArr[iTh], phiArr[iS], -1))

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

        self.t = self.t + self.delT
        
        # Compute an alpha shape 
        print("Computing alpha shape")
        self.ComputeHull()
        
        # Reconcile the shape to hull
        print("Reconciling alpha shape to hull")
        self.ReconcileHull()

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

        # Compute optimal control changes
        # See paper explanation (or Mathematica notebook 3D_Derivations.nb)
        if (phi < 0.001):
            # Avoid numerical blowup for degenerate phi
            thDot = 0
        else:
            thDot = -dUxdy*cos(th)**2 + dUydx*sin(th)**2            \
                +   (dUxdx-dUydy)*cos(th)*sin(th)                   \
                +   1.0 / tan(phi) * ( dUzdx*sin(th) - dUzdy*cos(th) )
        
        phiDot = -cos(phi)**2 * (dUzdx*cos(th) + dUzdy*sin(th)) \
                + sin(phi)**2 * (dUxdz*cos(th) + dUydz*sin(th)) \
                - 0.5 * ( -dUzdz + dUxdx*cos(th)**2 + dUydy*sin(th)**2 \
                        + (dUxdy + dUydx)*cos(th)*sin(th) ) * sin(2*phi)


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

        # TODO: Implement some kind of retry method to catch numerical
        # issues in the initial raycasting? Maybe try to cast rays that
        # aren't near any lines.

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

        p0 = np.array([p0x,p0y,p0z])
        p1 = np.array([p1x,p1y,p1z])

        # Find the triangular facet which intersects the line nearest p0
        firstTri = -1
        smallestR = -1
        for i in range(len(aShape)):

            # Represent the triangle as a list of points
            tri = aShape[i]
            triArr = np.array([pts[aShape[i,0],0:3], pts[aShape[i,1],0:3], pts[aShape[i,2],0:3]])

            # Compute the intersection parameter
            rInt = self.RayIntersectTriangle(p0, p1, triArr)
            if rInt == None:
                continue
           
            # Take the closest intersecter
            if firstTri == -1 or  rInt < smallestR:
                firstTri = i
                smallestR = rInt

        # The subset of the alpha shape which is in the hull
        hullFaces = set()
        # The facets to be processed as we traverse the list
        toProcess = Queue.Queue()

        # Re-order the points of this triangle (if needed) so that the
        # normal points outward
        ray = p1 - p0
        tNorm = self.TriNorm(aShape[firstTri])
        triNums = aShape[firstTri,:]
        if np.dot(ray, tNorm) > 0:
            swap = triNums[1]
            triNums[1] = triNums[2]
            triNums[2] = swap
        triNums = self.CanonTri(tuple(triNums))

        # Add this first triangle to the hull face set and the 
        # to-process list
        hullFaces.add(triNums)
        toProcess.put(triNums)

        # We will need a lookup table from points to facets to
        # efficiently implement the walk
        facetsForPoints = defaultdict(list)
        for i in range(len(aShape)):
            facetsForPoints[tuple(sorted((aShape[i,0],aShape[i,1])))].append(i)
            facetsForPoints[tuple(sorted((aShape[i,1],aShape[i,2])))].append(i)
            facetsForPoints[tuple(sorted((aShape[i,2],aShape[i,0])))].append(i)

        
        # Walk the facets to reconcile
        # Note: I haven't examined this closely enough to identify exactly which
        # kinds of topologically "wrong" alpha shapes will trick it in
        # to walking seemingly valid hull. See AlphaReconciliation.txt
        while not toProcess.empty():

            tri = toProcess.get()

            # The triangle was inserted in to toProcess in proper order,
            # so we can directly compute the normal vector
            currNVect = self.TriNorm(tri)

            # Check each of the neighbors
            for i in range(3):

                j = (i+1)%3
                k = (i+2)%3

                neighs = facetsForPoints[tuple(sorted((tri[i],tri[j])))]

                # If there's only one entry it's the current triangle,
                # and thus this facet has no neighbors.
                if len(neighs) == 1:
                    continue

                # Find the outermost neighbor
                outerMost = None
                outerMostTh = 999
                for oTriNum in neighs:

                    oTri = tuple(aShape[oTriNum])

                    # One of the elements in the list is the current
                    # facet, skip it
                    if self.CanonTri(oTri) == tri or self.CanonTri((oTri[0],oTri[2],oTri[1])) == tri:
                        continue

                    # Arrange the winding so the normal is consistent
                    # with the current triangle
                    for ioTri in range(3):
                        if oTri[ioTri] != tri[i] and oTri[ioTri] != tri[j]:
                            thirdPt = oTri[ioTri]

                    # The other triangle, properly wound
                    nOTri = (tri[i], thirdPt, tri[j])
                    # Take the opposite of the normal vector here so the
                    # math works out below
                    otherNVect = -1*self.TriNorm(nOTri)

                    # The vector along the edges, which defines the
                    # normal vector of the plane that we will take the angle
                    # about
                    pI = pts[tri[i],0:3]
                    pJ = pts[tri[j],0:3]
                    vEdge = pJ - pI

                    # Compute signed angle between the two normal
                    # vectors (one reversed) in the plane defined by vEdge
                    # This method explained here http://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane-reci
                    sinA = np.linalg.norm(np.cross(currNVect, otherNVect))
                    cosA = np.dot(currNVect, otherNVect)
                    theta = np.arctan2(sinA, cosA)
                    sign = np.dot(vEdge, np.cross(currNVect, otherNVect))
                    if sign < 0:
                        theta = -theta
                    if theta < 0:
                        theta += 2*np.pi

                    if theta < outerMostTh:
                        outerMost = nOTri
                        outerMostTh = theta

                # Now that the appropriate neighbor has been found,
                # rotate it until it looks like the initial
                # representation, with at most the last two switched
                cTri = self.CanonTri(outerMost)

                # Verify that the facet with the opposite normal vector
                # is not already there. If so, just move on because
                # this can happen legitmiately during joining
                # situations
                revTri = self.CanonTri((cTri[0], cTri[2], cTri[1]))
                if revTri in hullFaces:
                    continue

                # If this facet isn't already a hull facet, add it and
                # mark it for processing
                if cTri not in hullFaces:
                    # Note that this possibly adds the same triangle
                    # to hullFaces, ordered two different ways. This
                    # seems more proper anything else for now
                    hullFaces.add(cTri)
                    toProcess.put(cTri)
    
        # Finally, we have the set of all hull faces. Assign it and
        # return success
        print("Before reconciliation, alpha shape had %d faces, after it has %d facets"%(len(aShape),len(hullFaces)))
        self.surface = np.array(list(hullFaces))
            

    # Rotate the indicies defining a triangle to get a canonical
    # (comparable!) representation. Note that this DOES NOT change the
    # winding direction.
    def CanonTri(self, tri):
        while not (tri[0] < tri[1] and tri[0] < tri[2]):
            tri = (tri[2], tri[0], tri[1])
        return tri

    # Why doesn't numpy have normalize()?!
    def normalize(self,v):
        norm=np.linalg.norm(v)
        if norm==0: 
           return v
        return v/norm 

    # Takes triangle point indicies and computes the normal vector. Uses
    # the current pointset
    def TriNorm(self, triInds):
        pts = self.pointSets[-1]
        triArr = np.array([pts[triInds[0],0:3], pts[triInds[1],0:3], pts[triInds[2],0:3]])
        tNorm = np.cross(triArr[1] - triArr[0], triArr[2] - triArr[0])
        return self.normalize(tNorm)

    # Tests if a ray intersects a triangle. Returns none if there is no
    # intersection, or float s >= 0 which is a paramaterization of the
    # intersection point along p0 --> p1
    # 
    # p0 and p1 are points. tri should be a 3x3 array of point x coord
    #
    # From http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()
    # TODO do they have something I can cite? If they don't maybe use a
    # different algorithm with peer-review
    def RayIntersectTriangle(self, p0, p1, tri):
        
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

