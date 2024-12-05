from sympy import  symbols, cos, exp, I, pi, diff, lambdify, conjugate, simplify
import math

t, x, c, d = symbols('t x c d')

N=1000
h=1/1000
tau = 0.01
C=0.5
D=0.5

uTikslus =(t / 2 - exp(t) * I) * cos(4 * pi * x)

uTikslusFunc  = lambdify((x, t), uTikslus)

derivativeT = diff(uTikslus, t)
derivativeTFunc = lambdify((x, t), derivativeT)
print(derivativeTFunc(0,0))

derivativeX = diff(uTikslus, x)
derivativeXFunc = lambdify((x, t), derivativeX)
print(derivativeXFunc(0,1))

derivativeXSecond = diff(uTikslus, x, 2)
derivativeXSecondFunc = lambdify((x, t), derivativeXSecond)
print(derivativeXSecondFunc(0,0))

f = derivativeT - derivativeXSecond*I-I*c*uTikslus-I*d*(uTikslus*conjugate(uTikslus)*uTikslus)
print(f)


fsimplified = (-4*I*c*(t - 2*I*exp(t)) - I*d*(t - 2*I*exp(t))**2*(2*I*exp(conjugate(t)) + conjugate(t))*cos(4*pi*x)*cos(4*pi*conjugate(x)) + 64*I*pi**2*(t - 2*I*exp(t)) - 8*I*exp(t) + 4)*cos(4*pi*x)/8
fLambda = lambdify((x,t, c, d), fsimplified)

print(fLambda(0.5, 6, 0.5, 0.5))

res = (uTikslusFunc(1*h, tau+tau)-uTikslusFunc(1*h, tau))/tau
approx = 0.5j*((uTikslusFunc(2*h, tau+tau)-2*uTikslusFunc(1*h, tau+tau)+uTikslusFunc(0*h, tau+tau))/h*h +(uTikslusFunc(2*h, tau)-2*uTikslusFunc(1*h, tau)+uTikslusFunc(0*h, tau))/h*h)+C*1j*(uTikslusFunc(1*h, tau+tau)+ uTikslusFunc(1*h, tau))/2+D*1j*pow(abs((uTikslusFunc(1*h, tau+tau)+uTikslusFunc(1*h, tau))/2), 2)*((uTikslusFunc(1*h, tau+tau)+uTikslusFunc(1*h, tau))/2)+(fLambda(1*h, tau+tau, C, D)+fLambda(1*h, tau+tau, C, D))/2
print("--------------------------")
print(res)
print("--------------------------")
print(approx)