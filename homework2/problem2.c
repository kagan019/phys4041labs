//gcc C/tridag.c C/nrutil.c problem2.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//from Numerical Recipes Software(C) Copr. 1986-92. get it at http://numerical.recipes/
void tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n);

void gen_214tridiag(int n, float *ret[3]) {
    /* fills ret with 
        {{null,null,1,1,1,1,...},
        {null,2,4,4,4,...,2},
        {null,1,1,1,1,...,null}}
        each size n+2
    */
   for (int i = 0; i <= 2; i+=2) {
       for (int j = 1; j < n+2; j++) {
           ret[i][j] = 1;
       }
   }
   ret[1][1] = 2;
   ret[1][n+1] = 2;
   for (int j = 2; j < n+1;j++) {
       ret[1][j] = 4;
   } 
}

void gen_r(int n, float y[], float *ret) 
// n is the largest index of y (size-1)
{
    /*  fills ret with 
        {null,r0,r1,r2,...,rn)}
    */
   ret[1] = 3*(y[1] - y[0]);
   for (int i = 2; i <= n; i++) {
       ret[i] = 3*(y[i]-y[i-2]);
   }
   ret[n+1] = 3*(y[n] * y[n-1]);
}

void parametric_coeffs(float Di, float Dipo, 
    float yi, float yipo, float ret[4]) {
    /* ret is filled with the point-coefficients
        {ai,bi,ci,di}
    according to the wolfram page http://mathworld.wolfram.com/CubicSpline.html.
        a_i	=	y_i	
        b_i	=	D_i	
        c_i	=	3(y_(i+1)-y_i)-2D_i-D_(i+1)	
        d_i	=	2(y_i-y_(i+1))+D_i+D_(i+1)
    where y_i is {v[i],p[i]}
    */
    ret[0] = yi;
    ret[1] = Di;
    ret[2] = 3*(yipo-yi)-2*Di-Dipo;
    ret[3] = 2*(yi-yipo)+Di+Dipo;
}


void cubic_spline(int sz, float *v, float *p, float ret[][2][4]) 
//sz is the number of data points
{
    /* fills ret with
        {{x(t),y(t)}_0,{}_1,...,{..}_(sz-1)} parametric coefficients of spline by interval
    */
    float *bands[3];
    for (int i = 0; i < 3; i++) {
        bands[i] = malloc(sizeof(float)*(sz+1));
    }
    gen_214tridiag(sz-1,bands); //bands holding a,b,c for tridiag function
    float *vr = malloc(sizeof(float)*(sz+1));
    gen_r(sz-1,v,vr);
    float *pr = malloc(sizeof(float)*(sz+1));
    gen_r(sz-1,p,pr);
    float *vD = malloc(sizeof(float)*(sz+1));
    tridag(bands[0],bands[1],bands[2],vr,vD,sz);
    float *pD = malloc(sizeof(float)*(sz+1));
    tridag(bands[0],bands[1],bands[2],pr,pD,sz);
    for(int ivl = 0; ivl < sz-1; ivl++) {
        //ret[ivl] {xcoeffs, ycoeffs} along t of [0,1]
        parametric_coeffs(vD[ivl+1],vD[ivl+2],v[ivl],v[ivl+1],ret[ivl][0]); //xcoeffs
        parametric_coeffs(pD[ivl+1],pD[ivl+2],p[ivl],p[ivl+1],ret[ivl][1]); //ycoeffs 
    }
    for (int i = 0; i < 3; i++) {
        free(bands[i]);
    }
    free(vr);
    free(pr);
    free(vD);
    free(pD);
}

int interval_of_sign_change(float ptgt, float p[], int size) {
    int interval_of_sign_change = 0;
    for(;interval_of_sign_change < size-1; interval_of_sign_change++) {
        const int rising = p[interval_of_sign_change] < ptgt && ptgt < p[interval_of_sign_change+1];
        const int falling = p[interval_of_sign_change] > ptgt && ptgt > p[interval_of_sign_change+1];
        if (rising || falling) {
            break;
        }
    } 
    return interval_of_sign_change;
}

void offset(float ptgt, float p[], int sz) {
    for (int i = 0; i < sz; i++) {
        p[i] -= ptgt;
    }
    // so now the zeros of the spline are the solutions for y=ptgt
}

float eval_spline(float coef[4], float t) {
    return coef[0] + coef[1] * t + coef[2] * t*t
        + coef[3] * t * t * t;
}

int sign(float a) {
    return a<0;
}

void zero(float cubic_coeffs[4], float low, float high, float **ret) 
/*
    try to find a single zero. if it couldn't, set *ret to NULL
    [low, high] must be an 'interval_of_sign_change'
*/
{
    // bisection method
    //https://github.umn.edu/vinals/PHYS4041-2020/blob/0170af9eab30e8cab5b6a572e4abb95b51fdcaff/notes/function_zeros.ipynb  
    float a, b;
    a = low; b = high;
    const float eps = 0.00001;
    float ksi = (a+b)/2;
    for (int i = 0; i < 1000; i++) {
        if (fabs(eval_spline(cubic_coeffs,a)
            - eval_spline(cubic_coeffs,b)) < eps) 
            break;
        assert(sign(eval_spline(cubic_coeffs, a)) != sign(eval_spline(cubic_coeffs,b)));
        
        if (sign(eval_spline(cubic_coeffs,ksi)) == sign(eval_spline(cubic_coeffs,a))) {
            a = ksi;
        } else {
            b = ksi;
        }
        ksi = (a+b)/2.0;
    }
    if (fabs(eval_spline(cubic_coeffs,ksi)) > eps)
        *ret = NULL;
    else
        **ret = ksi;
}
    


int main() {
    const int size = 6;
    float v[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
    float p[] = {4.7, 3.5, 3.0, 2.7, 3.0, 2.4};
    const float ptgt = 3.25;

    //the cubic spline between the second and third points
    //is the only one that is guaranteed to cross 3.25.
    const int ioi = interval_of_sign_change(ptgt,p,size);
    offset(ptgt,p,size);
    // for this set of data, we expect 3.25 to be somewhere
    // between the second and third data point
    assert(ioi == 1); 

    float splines_by_interval[size-1][2][4]; // parametric coefficients for (x(t), y(t))
    cubic_spline(size,v,p,splines_by_interval);
    // thats part a. now b:
    float z; float *zz = &z;
    zero(splines_by_interval[ioi][1],0,1,&zz);
    if (zz == 0) {
        printf("no solutions :(\n");
        return 0;
    }
    printf("p=%f is interpolated to %6.6f\n", ptgt,
        eval_spline(splines_by_interval[ioi][0], z));
    return 0;
}
