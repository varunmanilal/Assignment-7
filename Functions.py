#List of functions
# 1.Gauss Jordan
# 2.Partial Pivot
# 3.Matrix multiplication
# 4.Forward substitution
# 5.Reverse substitution
# 6.LU decomposition
# 7.Reading Files
# 8.Reading 1D files
# 9.First derivative
# 10.Second derivative
# 11.Bracketing function
# 12.Bisection method
# 13.Regula Falsi method
# 14.Newton Raphson method
# 15.Laguerre's method
# 16.Synthetic division
# 17.Writing in notepad
# 18.Midpoint Method
# 19.Trapezoidal
# 20.Simpson
# 21.Monte Carlo
#22.Euler's method
#23.Runge Kutta 4
#24.Shooting method

def GaussJordan(A, B):
    # the pivot row
    pivot = A[k][k]
    for i in range(k, n):
        A[k][i] /= pivot
    B[k] = B[k] / pivot
    # other rows
    for i in range(n):
        if A[i][k] == 0 or i == k:
            continue
        else:
            term = A[i][k]
            for j in range(k, n):
                A[i][j] = A[i][j] - term * A[k][j]
            B[i] = B[i] - term * B[k]
    return B
#2
def Parpivot(A,B,k):
    if A[k][k] == 0:
        for i in range(k + 1, 4):
            if abs(A[i][k]) > abs(A[k][k]) and abs(A[k][k]) == 0:
                A[k], A[i] = A[i], A[k]
                B[k], B[i] = B[i], B[k]
    return A, B
#3
def MatrixMultiply(M,A):
    B=[]
    for i in range(len(M)):
        row =[]
        for j in range(len(A[0])):
            row.append(0)
        B.append(row)

    for x in range(len(M)):
        for y in range(len(A[0])):
            for z in range(len(M[0])):
                B[x][y] += M[x][z] * A[z][y]
    return B
#4
def fwrdsub(A , B):
    global Y
    Y = []
    for k in range(len(A)):
        Y.append(float(0))
    for i in range(0,len(A)):
        val = 0
        for j in range(0,i):
            val+=A[i][j]*Y[j]
        Y[i] += B[i] - val
    return Y
#5
def bkwdsub(A , B):
    global X
    X = []
    for k in range(len(A)):
        X.append(float(0))
    for i in reversed(range(len(A))):
        val = 0
        for j in reversed(range(0,len(A))):
            if j>i:
                val += A[i][j]*X[j]
        X[i] += (1/A[i][i])*(B[i] - val)
    return X
#6
def LUdecomp(A,C):
    for j in range(len(A)):
        Parpivot(A, C, j)
        for i in range(len(A)):
            if i<=j:
                sumt = 0
                for k in range(0, i):
                        sumt += A[i][k] * A[k][j]
                A[i][j] = A[i][j] - sumt
            if i>j:
                sumt =0
                for k in range(0, j):
                        sumt += A[i][k] * A[k][j]
                A[i][j] =(1/A[j][j])*(A[i][j]-sumt)

    print("The combined L-U matrix",A)
    return A
#7
def readfile(f1):
    a = f1.read()
    a1 = [item.split(' ') for item in a.split('\n')]
    A = []
    for i in range(len(a1)):
        row = []
        for j in range(len(a1[0])):
            row.append(0)
        A.append(row)
    for i in range(len(a1)):
        for k in range(len(a1[0])):
            A[i][k] = float(a1[i][k])
    return A
#8
def read1d(p):
    b = p.read()
    b1 = (b.split(' '))
    C=[]
    for i in range(len(b1)):
        C.append(0)
    for i in range(4):
        C[i] = float(b1[i])
    return C
#9
def firstderv(x,h,func):
    s = (func(x+h)-func(x-h))/(2*h)
    return s
#10
def secndderv(x,h,func):
    s = (func(x+h)+func(x-h)-2*func(x))/(h*h)
    return s


import math
#11
def bracketing(a,b,func):
    i=0
    #no of iterations limited to 12
    while func(a)*func(b)>0 and i<12:
        #if a or b is the root
        if func(a)*func(b)==0:
            if func(a)==0:
                print("the root is",a)
            if func(b)==0:
                print("The root is",b)
        elif func(a)*func(b)<0:
            continue
        #To move a and b
        else:
            if abs(func(a))<abs(func(b)):
                a-=1.5*(b-a)
            elif abs(func(a))>abs(func(b)):
                b+=1.5*(b-a)
        i+=1
    return a,b
#12
def bisection(a,b,func):
    #bracketing of the root
    (a,b)=bracketing(a,b,func)
    i=0
    #Y to store errors in every iteration
    Y = []
    print("Iterations", "  ","Absolute Error")
    while i<30 and abs(b-a)>0.000001:
        i+=1
        c = (a+b)/2
        if func(a)*func(c)<0:
            b=c
        elif func(a)*func(c)>0:
            a=c
        #if a or b turn out to be the root
        else:
            if func(a)==0:
                print("the root is",a)
            if func(c)==0:
                print("The root is",c)
        print(i,"           ",abs(b-a))
        Y.append(abs(b-a))
    return a,Y
#13
def regulafalsi(a,b,func):
    (a,b)=bracketing(a,b,func)
    i=0
    Y = []
    l=0
    c=b
    print("Iterations", "  ", "Absolute Error")
    while i<50 and abs(l-c)>0.000001:
        i+=1
        l=c
        #Regula falsi equation for c
        val =((b-a)*func(b))/(func(b)-func(a))
        c =b - val
        if func(a)*func(c)<0:
            b=c
        elif func(a)*func(c)>0:
            a=c
        #If either of them turn out to be the root
        else:
            if func(a) == 0:
                print("the root is", a)
            if func(c) == 0:
                print("The root is", c)
        print(i,"           ",abs(l-c))
        Y.append(abs(l-c))
    return c,Y
#14
def newtonraph(x,func):
    i = 0
    Y = []
    y=0
    print("Iterations", "  ", "Absolute Error")
    while i < 200 and abs(x-y) > 0.000001:
        i+=1
        y=x
        x = x - func(x)/firstderv(x,0.01,func)
        Y.append(abs(x-y))
        print(i,"           ",abs(x-y))
    return x,Y
#15
def laguerre(alpha,n,func):
    Y = []
    if func(alpha)!=0:
        alpha1=alpha+1
        i=0
        while i<50 and abs(alpha1-alpha)>0.000001:
            i+=1
            #to find difference between consecutive iterations
            alpha1=alpha
            #value of G and H
            G = firstderv(alpha,0.1,func)/func(alpha)
            H = G*G-(secndderv(alpha,0.1,func)/func(alpha))

            if abs(G+math.sqrt((n-1)*(n*H-math.pow(G,2))))<abs(G-math.sqrt((n-1)*(n*H-math.pow(G,2)))):
                a = n/(G-math.sqrt((n-1)*(n*H-math.pow(G,2))))
            else:
                a = n/(G+math.sqrt((n-1)*(n*H-math.pow(G,2))))
            alpha = alpha -a
            Y.append(abs(alpha1 - alpha))
        print(round(alpha,6),"is one of the roots")
    else:
        #If alpha called in main function is one of the roots
        print(round(alpha,6),"is one of the roots")

    return alpha,Y
#16
def syndiv(coeffs,alpha):
    #Making absolute value of 1st coefficient 1
    if abs(coeffs[0])!=1:
        for k in range(len(coeffs)):
            coeffs[k]=coeffs[k]/coeffs[0]

    for l in range(1,len(coeffs)):
        coeffs[l]+=alpha*coeffs[l-1]
    print("The coefficient matrix is ",coeffs)
    return coeffs
#17
#To print in notepad
def Notepadwrt(q1,y1):
    j = 1
    q1.write('Itr       Value \n')
    for i in y1:
        q1.write(str(j) + '         ')
        q1.write(str(i) + '\n')
        j += 1
    q1.close()

#The midpoint function
def midpoint(a,b,N,func):
    #The interval value
    h = (b-a)/N
    sum=0
    for i in range(N):
        value = (func(a+i*h)+func(a+(i+1)*h))/2
        sum+=h*value
    return sum

#The trapezoidal function
def trapezoidal(a,b,N,func):
    h = (b-a)/N
    sum=0
    #There are 2 weigth functions
    for i in range(N+1):
        if i==0 or i==N:
            sum+=h*func(a+i*h)/2
        else:
            sum+=h*func(a+i*h)
    return sum

#Simpson's method
def simpson(a,b,N,func):
    h = (b - a) / N
    sum = 0
    #There are 3 weigth functions
    for i in range(N + 1):
        if i == 0 or i == N:
            sum += h * func(a + i * h) / 3
        elif i%2 == 1:
            sum += 4 * h * func(a + i * h) / 3
        else:
            sum += 2 * h * func(a + i * h) / 3
    return sum

#Monte carlo method
import random
#The limits ,no of iterations and the func is the input
def montecarlo(a,b,m,func):
    j = 0
    total =0
    #Using random function
    while j<m:
        x = random.uniform(a, b)
        ans = (b-a)/m*func(x)
        total+=ans
        j+=1
    return total
#Euler method with inputs,x0 , y0, range, interval and function
def Eulerexp(x,y,n,h,f):
    #print("Values obtained for h =", h)
    #Values stored in arrays
    X=[]
    Y=[]
    #calculate when x less than range
    while x<=n:
        y = y + h*f(x,y)
        x = x + h
        Y.append(y)
        X.append(x)
    return X,Y
#Runge Kutta4 with inputs x0,y0,z0 ,interval, range for x and both the functios
def RK4(x,y,v,h,r,funcf,funcg):
    #Values stored in arrays
    X=[]
    Y=[]
    #Calculate when x less than range
    while x <= r:
        k1 = h*funcf(v)
        l1 = h*funcg(x,y,v)

        k2 = h * funcf(v+l1/2)
        l2 = h * funcg(x+h/2,y+k1/2, v+l1/2)

        k3 = h * funcf(v + l2 / 2)
        l3 = h * funcg(x + h / 2, y + k2 / 2, v + l2 / 2)

        k4 = h * funcf(v + l3)
        l4 = h * funcg(x + h, y + k3, v + l3)
        #Calculate for y and v i.e z for next iteration
        y = y + 1/6*(k1+2*k2+2*k3+k4)
        v = v + 1/6*(l1+2*l2+2*l3+l4)
        x = x + h
        X.append(x)
        Y.append(y)
    return X,Y
#RK4 + Runge kutta
def shootingmethod(funcf,funcg,yn):
    #X and Y for 2 inputs
    def RungeKutta(zh,zt):
        (X1, Y1) = RK4(0, 1, zh, 0.01, 1, funcf, funcg)
        (X2, Y2) = RK4(0, 1, zt, 0.01, 1, funcf, funcg)
        return X1,Y1,X2,Y2
    #Inputs by user for 2 values
    def estimate():
        for i in range(20):
            est1 = input("Estimate 1st value for slope ")
            zh = float(est1)
            est2 = input("Estimate 2nd value for slope ")
            zt = float(est2)
            (X1, Y1, X2, Y2) = RungeKutta(zh, zt)
            #If both values are on same side of yn then take a diff input
            if (Y2[99]>yn and Y1[99]>yn) or (Y2[99]<yn and Y1[99]<yn):
                print("Try another set")
                continue
            else:
                return zh, zt
    (zh,zt) = estimate()
    for i in range(20):
        (X1, Y1, X2, Y2) = RungeKutta(zh, zt)
        #If values are close enough to yn stop loop
        if abs(Y1[99]-yn)<0.002:
            print("The solution for z is :",zh)
            return Y1[99],Y1,X1
        if abs(Y2[99]-yn)<0.002:
            print("The solution for z is :",zt)
            return Y2[99],Y2,X2
        #accounting for both input cases, y1<Y2 and y1<y2
        if Y1[99]>yn and Y2[99]<yn:
            z=zt+(((zh-zt)/((Y1[99])-Y2[99]))*((yn)-Y2[99]))
            Y3,X3=RK4(0, 1, z, 0.01, 1, funcf, funcg)
            if abs(Y3[99]-yn)<0.02:
                print("The obtained solution for z is :",z)
                return Y3[99],Y3,X3
            elif z<yn:
                zt=z
            else:
                zh=z
        if Y2[99]>yn and Y1[99]<yn:
            z=zh+(((zt-zh)/((Y2[99])-Y1[99]))*((yn)-Y1[99]))
            Y3,X3=RK4(0, 1, z, 0.01, 1, funcf, funcg)
            if abs(Y3[99]-yn)<0.02:
                print("The obtained solution for z is :",z)
                return Y3[99],Y3,X3
            elif z<yn:
                zh=z
            else:
                zt=z