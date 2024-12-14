from sympy import sin,  symbols, cos, exp, I, pi, diff, lambdify, conjugate, simplify
import math

t, x, c, d = symbols('t x c d')

N=1000
h=1/N
tau = 0.01
C=0.5
D=0.5

uTikslus =(t / 2 - exp(t) * I) * cos(4 * pi * x)

uTikslusFunc  = lambdify((x, t), uTikslus)

derivativeT = diff(uTikslus, t)
derivativeTFunc = lambdify((x, t), derivativeT)
#print(derivativeTFunc(0,0))

derivativeX = diff(uTikslus, x)
derivativeXFunc = lambdify((x, t), derivativeX)
#print(derivativeXFunc(0,1))

derivativeXSecond = diff(uTikslus, x, 2)
derivativeXSecondFunc = lambdify((x, t), derivativeXSecond)
#print(derivativeXSecondFunc(0,0))

f = derivativeT - derivativeXSecond*I-I*c*uTikslus-I*d*(uTikslus*conjugate(uTikslus)*uTikslus)
print(f)

#print(derivativeX)
#print(derivativeXFunc(0,1))
#print(derivativeXFunc(1,1))
#print(sin(4 * pi).evalf())
#print(derivativeX.subs({x: 1, t: 1}))


#fsimplified = (-4*I*c*(t - 2*I*exp(t)) - I*d*(t - 2*I*exp(t))**2*(2*I*exp(conjugate(t)) + conjugate(t))*cos(4*pi*x)*cos(4*pi*conjugate(x)) + 64*I*pi**2*(t - 2*I*exp(t)) - 8*I*exp(t) + 4)*cos(4*pi*x)/8
fLambda = lambdify((x,t, c, d), f )

#print(fLambda(0.5, 6, 0.5, 0.5))

res = (uTikslusFunc(1*h, tau+tau)-uTikslusFunc(1*h, tau))/tau
def uj(x, t):
    return uTikslusFunc(x, t+tau)
def ujplus1(x, t):
    return uTikslusFunc(x+h, t+tau)
def ujminus1(x, t):
    return uTikslusFunc(x-h, t+tau)

def u(x, t):
    return uTikslusFunc(x, t)
def uplus1(x, t):
    return uTikslusFunc(x+h, t)
def uminus1(x, t):
    return uTikslusFunc(x-h, t)

def algorithmPart1(x, t):
    return 0.5j*((ujplus1(x, t)-2*uj(x, t)+ujminus1(x, t))/(h*h) +(uplus1(x, t)-2*u(x, t)+uminus1(x, t))/(h*h))

def algorithmPart2(x, t):
    return C*1j*(uj(x, t)+ u(x, t))/2+D*1j*pow(abs((uj(x, t)+ u(x, t))/2), 2)*((uj(x, t)+ u(x, t))/2)

def algorithmPart3(x, t):
    return (fLambda(x, t+tau, C, D)+fLambda(x, t, C, D))/2

def algorithm(x, t):
    approx = algorithmPart1(x, t)+algorithmPart2(x, t)+algorithmPart3(x, t)
    return approx

def simple(x, t):
    return (uj(x, t)-u(x, t))/tau
#approx = 0.5j*(()/h*h +(uTikslusFunc(2*h, tau)-2*uTikslusFunc(1*h, tau)+uTikslusFunc(0*h, tau))/h*h)+C*1j*(uTikslusFunc(1*h, tau+tau)+ uTikslusFunc(1*h, tau))/2+D*1j*pow(abs((uTikslusFunc(1*h, tau+tau)+uTikslusFunc(1*h, tau))/2), 2)*((uTikslusFunc(1*h, tau+tau)+uTikslusFunc(1*h, tau))/2)+(fLambda(1*h, tau+tau, C, D)+fLambda(1*h, tau, C, D))/2


#print(simple(0.31, 1.7)-algorithm(0.31, 1.7))
#print(simple(0.031, 0.17)-algorithm(0.031, 0.17))
#print((simple(0.31, 1.7)-algorithm(0.31, 1.7))/(simple(0.031, 0.17)-algorithm(0.031, 0.17)))
print("----------test1----------")
res1=abs(simple(0.1, 1.8)-algorithm(0.1, 1.8))
h=h/10 
tau=tau/10
res2=abs(simple(0.1, 1.8)-algorithm(0.1, 1.8))
print(res1/res2)




############## test2
print("----------test2----------")
h=h*10
tau=tau*10

def Cconst(h, tau):
 return (2-(2j*pow(h, 2))/tau)

def cPart(x, t, h, tau):
    return ujplus1(x, t)-uj(x, t)*Cconst(h, tau)+ujminus1(x, t)

def fPart(x, t):
    return (uplus1(x, t)-2*u(x, t)+uminus1(x, t)-1j*(2*h*h)/tau*u(x, t)-2j*h*h*algorithmPart2(x, t)-1j*h*h*(fLambda(x, t+tau, C, D)+fLambda(x, t, C, D)))




restest2prev=abs(cPart(0.71, 3, h, tau)+fPart(0.71, 3))
h=h/10
tau=tau/10
restest2=abs(cPart(0.71, 3, h, tau)+fPart(0.71, 3))
print(restest2prev/restest2)
print(restest2prev)
print(restest2)