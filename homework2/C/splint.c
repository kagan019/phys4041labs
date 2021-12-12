# include "nrutil.h"


void splint(float xa[],float ya[], float y2a[], int n, float x, float *y)
/*

Given arrays xa[1..n], ya[1..n] which tabulate the function (with the
xa's in order) and given the array y2a[1..n] which is the output from 
the function spline, and given a value of x, this routine returns the
cubic-spline interpolated value y

*/

{
void nrerror(char error_text[]);
int klo,khi,k;
float h,b,a;

klo=1;
khi=n;
while(khi-klo > 1) {
	k=(khi+klo) >> 1;
	if(xa[k] > x) khi=k;
	else klo=k;
}
h=xa[khi]-xa[klo];
if (h == 0) nrerror("Bad xa input to routine splint");
a=(xa[khi]-x)/h;
b=(x-xa[klo])/h;
*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
