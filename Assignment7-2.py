from Functions import *
import math
#1st coupled eqn
def funcf(z):
    p = z
    return p
#2nd coupled eqn
def funcg(x,y,z):
    p = 1-x-z
    return p
#calculation by calling RK4
(X,Y)=RK4(0,2,1,0.1,5,funcf,funcg)
#Printing of values for printing
for n in range(len(X)):
    print(X[n],Y[n])

#No results to be appended
