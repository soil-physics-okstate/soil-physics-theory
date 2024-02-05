#PSP_simpson.py
from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from math import *

def plottrapzd(func,a,b,n):
	it = 1; 
	for j in range(1,n):
		it <<= 1;
	x = np.zeros(2+2*it);
	y = np.zeros(2+2*it);
	d = (b-a)/float(it);
	x[0]=a;
	y[0]=func(a);
	x[1]=a+0.5*d;
	y[1]=func(a);
	for i in range(1,it):
		x[2*i]=a+(float(i)-0.5)*d;
		y[2*i]=func(a+float(i)*d);
		x[2*i+1]=a+(float(i)+0.5)*d;
		y[2*i+1]=func(a+float(i)*d);
	x[2*it]=b-0.5*d;
	y[2*it]=func(b);
	x[2*it+1]=b;
	y[2*it+1]=func(b);
	x_p = np.zeros(1+it);
	y_p = np.zeros(1+it); 
	for i in range(0,1+it):
		x_p[i] = x[2*i]+0.5*d;
		y_p[i] = y[2*i];
	x_p[0] = x[0];
	plt.fill_between(x,y,0,color='0.5')
	x2 = np.zeros(1000);
	y2 = np.zeros(1000);
	for i in range(0,1000):
		x2[i]=a+(b-a)/1000.*float(i+1);
		y2[i]=func(x2[i]);
	plt.plot(x2,y2,color='black')
	plt.ylim(0,max(y2)+1.);
	plt.xlim(a,b);
	plt.xlabel('x',fontsize=16,labelpad=4)
	plt.ylabel('y',fontsize=16,labelpad=4)
	plt.tick_params(axis='both', which='major', labelsize=16,pad=4)
	plt.tick_params(axis='both', which='minor', labelsize=16,pad=4)

def trapzd(func, a, b, n):
	if (n == 1): 
		trapzd.s=0.5*(b-a)*(func(a)+func(b));
		return trapzd.s;
	else:
		it = 1; 
		for j in range(1,n-1):
			it <<= 1;
		tnm=float(it);
		del_=(b-a)/tnm;
		x=a+0.5*del_;
		sum_ = 0.0;
		for j in range(1,it+1): 
			sum_ += func(x);
			x+=del_;
		trapzd.s=0.5*(trapzd.s+(b-a)*sum_/tnm);
		return trapzd.s;

def qsimp(func, a, b):
	EPS = 1.0e-6;
	JMAX = 20;
	ost=0.0;
	os=0.0;
	for j in range(1,JMAX+1):
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (j > 5):
			if (fabs(s-os) < EPS*fabs(os) or (s == 0.0 and os == 0.0)): 
				plottrapzd(func,a,b,j);
				plt.show() 
				return s;
		os=s;
		ost=st;
	return 10.0;

def func(x):
	#return (x-1.)*(x-1.)+0.6*x*sin(6.*x)+1.4;
	return (1/x);

def main():
	integral = qsimp(func,1,4);
	print ("Numerical solution:")
	print (integral)
main()
