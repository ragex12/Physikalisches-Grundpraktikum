from numpy import loadtxt, pi, linspace, sqrt, arange
import matplotlib.pylab as plt
from functions import linreg

def aufgabe1():
	f = loadtxt("../data/rl-kreis.dat")[:,0]
	I = loadtxt("../data/rl-kreis.dat")[:,1]
	U = loadtxt("../data/rl-kreis.dat")[:,2]
	phi = loadtxt("../data/rl-kreis.dat")[:,3]
	sigma_U = 3.0
	sigma_I = 2.0
	sigma_Z_2 = sqrt(sigma_U**2*4*U**2/I**4-sigma_I**2*4*U**4/I**6)
	Z_2 = U**2/I**2
	print(Z_2)
	omega_2 = ((2*pi*f)**2)*10**(-5)
	plt.figure(1)
	plt.xlabel(r"$\omega^2$[$10^{-5}\rm{Hz}^2$]")
	plt.ylabel(r"$Z^2$[$\Omega^2$]")
	[b1,m1,sigma_b1,sigma_m1] = linreg(omega_2,Z_2,sigma_Z_2)
	plt.errorbar(omega_2,Z_2,sigma_Z_2,fmt='ko')
	x = arange(0,90,1)
	plt.plot(m1*x+b1)
	plt.savefig("../latex/aufgabe1.pdf")
	
def aufgabe2():
	f = loadtxt("../data/rlc-kreis.dat")

