from numpy import loadtxt, pi, linspace, sqrt, arange
import matplotlib.pylab as plt
from functions import linreg, close2zero

sigma_U = 4.3
sigma_I = 3.3*10**-3

def aufgabe1():
	f = loadtxt("../data/rl-kreis.dat")[:,0]
	I = loadtxt("../data/rl-kreis.dat")[:,1]*10**-3
	U = loadtxt("../data/rl-kreis.dat")[:,2]
	sigma_Z_2 = sqrt(sigma_U**2*4*U**2/I**4-sigma_I**2*4*U**4/I**6)
	Z_2 = U**2/I**2
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
	f = loadtxt("../data/rlc-kreis.dat")[:,0]
	I = loadtxt("../data/rlc-kreis.dat")[:,1]*10**-3
	U_ges = loadtxt("../data/rlc-kreis.dat")[:,2]
	Z = U_ges/I
	omega = 2*pi*f
	sigma_Z = sqrt(sigma_U**2*(1/I**2)-sigma_I**2*(U_ges/I**2)**2)
	R_S = min(Z)
	for i in range(0,len(Z)):
		if(Z[i] == R_S):
			omega_R = omega[i]
	print([R_S, omega_R])
	plt.figure(2)
	plt.xlabel(r"$\omega$[Hz]")
	plt.ylabel(r"$Z$[$\Omega$]")
	plt.errorbar(omega,Z,sigma_Z,fmt='ks')
	plt.savefig("../latex/aufgabe2.pdf")


def aufgabe3():
	phi = loadtxt("../data/rlc-kreis.dat")[:,4]
	f = loadtxt("../data/rlc-kreis.dat")[:,0]
	omega = 2*pi*f
	sigma_phi = 2
	[omega_R, phi_0]= close2zero(omega,phi)
	print([omega_R, phi_0])
	plt.figure(3)
	plt.xlabel(r"$\omega$[Hz]")
	plt.ylabel(r"$\varphi$[°]")
	plt.errorbar(omega, phi, sigma_phi, fmt='kx')
	plt.savefig("../latex/aufgabe3.pdf")

aufgabe1()
aufgabe2()
aufgabe3()
