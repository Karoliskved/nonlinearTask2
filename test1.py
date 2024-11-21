import math

a=0
c=0.5
d=0.5
k=0
alpha=0
beta=0
gamma=0
theta=0


def u(x, t):
    return (t/2-math.exp(t)*1j)*math.cos(4*math.pi*x)

def derivativeT(x, t):
    return (0.5-(1+(math.pi*1j)/2)*math.exp(t)*1j)*math.cos(4*math.pi*x)

def derivativeX(x, t):
    return 2*math.pi(-t+2*math.exp(t)*1j)*math.sin(4*math.pi*x)

def secondDerivative(x, t):
    return 8*pow(math.pi, 2)*(-t+2*math.exp(t)*1j)*math.cos(4*math.pi*x)

def f(x, t):
    return derivativeT(x, t)-1j*secondDerivative(x, t)-1j*c*u(x, t)-1j*d*pow(abs(u(x,t)), 2)*u(x, t)

print(f(0.5, 5))