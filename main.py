import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# 
# Fourier Series Code
# 

#  funtions for fourier approximation
def f1(x):
    ''' f1(x) = x
    ''' 
    return x

def f2(x):
    ''' Square wave 
               ~ -1, -1 < x < 0
        f2(x)= |
               ~  1,  0 < x < 1
    '''
    return signal.square(2*x)

def f3(x):
    ''' Cubic
    '''
    
    return 5*x*(x+1.)*(x-1.)

def test(x):
    return x*(1-x)
    # return 2*x/x

def test2(x):
    return np.sin(2*np.pi*x)

def test3(x):
    return f1(x)+test2(x)+f3(x)

def f4(x):
    return 2*np.ones(len(x))

def test4(x):
    return np.sin(2*np.pi*x)+1.1*np.sin(4*np.pi*x)+np.sin(x)+np.sin(3*x)

def genSln(n, bn, L, x, t):
    return bn*np.exp(-(n*np.pi*alpha)**2*t)*np.sin(n*np.pi*x/L)

def fourierCoefs(x, nMax, x0, x1, f, A, B):
    """ Produces Fourier Series approximation of a function
    
        ----------

        Parameters

        ----------

        x: array_like
            sample values

        nMax: int
            number of terms

        x0: float
            lower bound on window/approximation

        x1: float
            upper bound on window/approximation

        f: function
            function to be approximated

    """
    T = x1 - x0

    # need new range over the actual range of the approximation (not the plotting range)
    # so integration is correct
    xFuncPeriod = np.arange(x0, x1, 0.01)  

    b = []

    # (f(xFuncPeriod)-A*xFuncPeriod-B)

    for n in range(nMax):
        # an = (2/T)*np.trapz(f(xFuncPeriod)()*np.cos(2*n*np.pi*xFuncPeriod*(1/T)), x=xFuncPeriod)
        bn = (2/T)*np.trapz((f(xFuncPeriod)-A*xFuncPeriod-B)*np.sin(2*n*np.pi*xFuncPeriod*(1/T)), x=xFuncPeriod)

        b.append(bn)

    return b

def update(frame, x, uFrames):
    ln.set_data(x,uFrames[:,int(frame)])
    return ln,

def steadySln(bc, L):
    B=bc[0]
    A=-(bc[0]-bc[1])/L

    return A,B

if __name__=="__main__":

    alpha=1.3e-2

    L=1
    xPoints=40
    tMax=1000

    pps=100
    tPoints=int(pps*tMax)

    bc=[0,4]

    x=np.linspace(0,L,xPoints)
    T=np.linspace(0,tMax,tPoints)

    A, B = steadySln(bc, L)
    uSteady=A*x+B

    b = fourierCoefs(x, 20, 0, L, test4, A, B)


    uFrames=np.zeros((xPoints, tPoints))

    for j, t in enumerate(T):
        u=np.zeros(len(x))

        for n, an in enumerate(b):
            u+=genSln(n,an,L,x,t)

        uFrames[:,j]=u+uSteady

    fig, ax = plt.subplots()
    ln, = plt.plot([], [])
    ax.set_xlim(0,L)
    ax.set_ylim(0,5)

    ani = FuncAnimation(fig, update, interval=20, frames=np.linspace(0, tPoints-1), fargs=(x, uFrames), blit=True)
    plt.show()


    