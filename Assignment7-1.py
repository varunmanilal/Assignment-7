from Functions import *
import math
#1st diff equation
def func(x,y):
    p = y*math.log(y,2.71828)/x
    return p
#2nd diff equation
def func2(x,y):
    p = 6 - (2*y)/x
    return p
#Different values of h for euler's method
H=[0.01,0.05,0.1,0.5]

for h in H:
    #Taken upto range 10
    (x,y)=Eulerexp(2,2.71,10,h,func)

for h in H:
    #taken upto range 10
    Eulerexp(3,1,10,h,func2)
#Values printed for graphs
#No results to be appended