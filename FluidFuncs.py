# Nicholas Sharp - nsharp3@vt.edu

# Various 2D fluid functions and derivatives for dynamical systems problems

import numpy as np

#  A a linear flow in the postive X direction
class LinearXFlow:
    
    def __init__(self, vel):
        self.vel = vel

    # Velocity in the x-direction
    def U(self,x,y,t):
        return self.vel

    # Velocity in the y-direction        
    def V(self,x,y,t):
        return 0

    # dU/dx   
    def Ux(self,x,y,t):        
        return 0
    
    # dU/dy   
    def Uy(self,x,y,t):        
        return 0

    # dV/dx   
    def Vx(self,x,y,t):        
        return 0

    # dV/dy   
    def Vy(self,x,y,t):        
        return 0

#  A 3D linear flow in the postive X direction
class LinearXFlow3D:
    
    def __init__(self, vel):
        self.name = 'LinearXFlow3D'
        self.info = '%s: vel = %.4e'%(self.name,vel)
        self.vel = vel

    # Velocity in the x-direction
    def Ux(self,x,y,z,t):
        return self.vel

    # Velocity in the y-direction        
    def Uy(self,x,y,z,t):
        return 0

    # Velocity in the z-direction        
    def Uz(self,x,y,z,t):
        return 0
    
    def dUxdx(self,x,y,z,t):        
        return 0
    
    def dUxdy(self,x,y,z,t):        
        return 0

    def dUxdz(self,x,y,z,t):        
        return 0

    def dUydx(self,x,y,z,t):        
        return 0
    
    def dUydy(self,x,y,z,t):        
        return 0

    def dUydz(self,x,y,z,t):        
        return 0

    def dUzdx(self,x,y,z,t):        
        return 0
    
    def dUzdy(self,x,y,z,t):        
        return 0

    def dUzdz(self,x,y,z,t):        
        return 0

class DoubleVortex:

    def __init__(self, vel, b, omega):
        self.name = 'DoubleVortex'
        self.info = '%s: vel = %.4e b = %.4e omega = %.4e'%(self.name,vel,b,omega)
        self.vel = vel
        self.b = b
        self.omega = omega

    def Ux(self,x,y,t):
    	arg = np.pi * (x + self.b * np.cos(self.omega * t))
    	return self.vel*np.sin(arg)*np.cos(np.pi*y)

    def Uy(self,x,y,t):
	    arg = np.pi * (x + self.b * np.cos(self.omega * t))
	    return -self.vel*np.cos(arg)*np.sin(np.pi*y)

    def dUxdx(self,x,y,t):
    	arg = np.pi * (x + self.b * np.cos(self.omega * t))
    	return self.vel*np.pi*np.cos(arg)*np.cos(np.pi*y)

    def dUydy(self,x,y,t):
    	arg = np.pi * (x + self.b * np.cos(self.omega * t))
    	return -self.vel*np.pi*np.cos(arg)*np.cos(np.pi*y)
    
    def dUxdy(self,x,y,t):
    	arg = np.pi * (x + self.b * np.cos(self.omega * t))
    	return -self.vel*np.pi*np.sin(arg)*np.sin(np.pi*y)

    def dUydx(self,x,y,t):
    	arg = np.pi * (x + self.b * np.cos(self.omega * t))
    	return self.vel*np.pi*np.sin(arg)*np.sin(np.pi*y)


class Gyre3D:

    def __init__(self, vel, a, b, c, omega):
        self.name = 'Gyre3D'
        self.info = '%s: vel = %.4e a = %.4e b = %.4e c = %.4e omega = %.4e'%(self.name,vel,a,b,c,omega)

        self.vel = vel
        self.a = a
        self.b = b
        self.c = c
        self.omega = omega

    # Velocity in the x-direction
    def Ux(self,x,y,z,t):
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return -self.vel * np.sin(valX) * np.cos(valY) * np.sin(valZ)

    # Velocity in the y-direction        
    def Uy(self,x,y,z,t):
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return -self.vel * np.cos(valX) * np.sin(valY) * np.sin(valZ)

    # Velocity in the z-direction        
    def Uz(self,x,y,z,t):
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return self.vel * np.cos(valX) * np.sin(valY) * np.cos(valZ)
    
    def dUxdx(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return - self.a * self.vel * np.cos(valX) * np.cos(valY) * np.sin(valZ)
    
    def dUxdy(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return self.b * self.vel * np.sin(valX) * np.sin(valY) * np.sin(valZ)

    def dUxdz(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return - self.c * self.vel * np.sin(valX) * np.cos(valY) * np.cos(valZ)

    def dUydx(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return self.a * self.vel * np.sin(valX) * np.sin(valY) * np.sin(valZ)
    
    def dUydy(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return -self.b * self.vel * np.cos(valX) * np.cos(valY) * np.sin(valZ)

    def dUydz(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return - self.c * self.vel * np.cos(valX) * np.sin(valY) * np.cos(valZ)

    def dUzdx(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return - self.a * self.vel * np.sin(valX) * np.sin(valY) * np.cos(valZ)
    
    def dUzdy(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return self.b * self.vel * np.cos(valX) * np.cos(valY) * np.cos(valZ)

    def dUzdz(self,x,y,z,t):        
        valX = self.a*x + self.omega*t
        valY = self.b*y + self.omega*t
        valZ = self.c*z + self.omega*t
        
        return - self.c * self.vel * np.cos(valX) * np.sin(valY) * np.sin(valZ)
