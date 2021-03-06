#!/usr/bin/env python3

from numpy import array, sqrt, arange, average

def weighted_mean(arr,err):
	#print(arr,err)
	a,b = average(arr,weights=1./err**2,returned=True)
	sigma_a = sqrt(1/b)
	return a,sigma_a

def linreg(x,y,yerr):
	if type(yerr) == float:
		cache=[]
		for i in arange(0,len(y)):
			cache.append(yerr)
		yerr = array(cache)
	
	sigma_sum = sum(1/yerr**2)
	x2_sigma_sum = sum(x**2/yerr**2)
	x_sigma_sum = sum(x/yerr**2)
	xy_sigma_sum = sum(x*y/yerr**2)
	y_sigma_sum = sum(y/yerr**2)
	Delta = sigma_sum*x2_sigma_sum-x_sigma_sum**2
	m=(1/Delta)*((sigma_sum*xy_sigma_sum)-(x_sigma_sum*y_sigma_sum))
	b=(1/Delta)*((x2_sigma_sum*y_sigma_sum)-(x_sigma_sum*xy_sigma_sum))
	sigma_m = sqrt((1/Delta)*sigma_sum)
	sigma_b = sqrt((1/Delta)*x2_sigma_sum)
	return [b,m,sigma_b,sigma_m]

def convert(data, sigma_data):
	str_data = format(data, '.10f').rstrip('0')
	str_sigma_data = format(sigma_data,'.10f').rstrip('0')
	a,b = str_sigma_data.split(".")
	s = (a+b).split("0")
	zeroing = 0
	for i in s:
		if i=='':
			zeroing = zeroing + 1
		else:
			break
	a,b = str_data.split(".")
	if len(b) < zeroing:
		for i in arange(zeroing - len(b)):
			b = b + '0'
	if str_sigma_data[zeroing+1] == '9':
		s = a+'.'+b[:zeroing-1]+'('+str(int(str_sigma_data[zeroing+1])-8)+')'
	else:
		s = a+'.'+b[:zeroing]+'('+str(int(str_sigma_data[zeroing+1])+1)+')'
	return s

#looking for P(x_0,y_0) closest to zero
def close2zero(x,y):
	y_0 = min(abs(y))
	x_0 = []
	for i in range(0,len(y)):
		if(y[i]==y_0):
			x_0.append(x[i])
	return [x_0,y_0]

