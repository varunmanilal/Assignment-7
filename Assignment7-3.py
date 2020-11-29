from Functions import *
import math
#1st coupled
def funcf(z):
    p = z
    return p

def funcg(x,y,z):
    p = z+1
    return p
#Value of yn given
yn = 3.43656

#Values stored for graph
a,b,c=shootingmethod(funcf,funcg,3.43)

#Result appended
#Estimate 1st value for slope 0.5
#Estimate 2nd value for slope 1.5
#The solution for z is : 0.9961801048227649






