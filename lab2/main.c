#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 500.
#define beta 0.001
#define lambda 0.1
#define TOL 1e-6
#define a11  0.25
#define a12 (0.25 - sqrt(3) / 6.)
#define a21 (0.25 + sqrt(3) / 6.)
#define a22  0.25
double fun(double u)
{
    double alfa = beta * N - lambda; 
    return alfa * u - beta * u * u;
}
double dfun(double u)
{
    double alfa = beta * N - lambda; 
    return alfa - 2. * beta * u;
}
void MetodaPicarda(FILE * fp)
{
    double tmax = 100.;
    double deltaT = 0.1;
    double u0 = 1.;
    double u = u0;
    double uq = 0.;
    double un = u0;
    int q  = 0;
    int qMax = 0; 
    for(double t = deltaT; t <= tmax; t+=deltaT)
    {     
        uq = 0.;
        while(fabs(u - uq) > TOL && q <= 20)
        {
            uq = u; 
            u = un + (deltaT/2.) * (fun(un) + fun(uq));      
            q++;         
        }
        if(q > qMax) 
            qMax = q;          
        un = u;   
        fprintf(fp,"%g %g %g %d\n",t,un, N - un,q); 
        q = 0;
    }
    printf("%d\n",qMax);
}
void MetodaNewtona(FILE * fp)
{
    double tmax = 100;
    double deltaT = 0.1;
    double u = 1.;
    double uq = 0.;
    double un = 1.;
    int q  = 0;
    int qMax = 0; 
    for(double t = deltaT; t <= tmax; t += deltaT)
    {
        uq = 0.;
        while(fabs(u - uq) > TOL && q <= 20)
        {
            uq = u;
            u = uq - (uq - un - (deltaT/2.) * (fun(un) + fun(uq))) / (1 - (deltaT/2.) * dfun(uq));
            q++;
        }
        if(q > qMax) 
            qMax = q;
        un = u;
        q = 0;
        fprintf(fp,"%g %g %g\n",t,un, N - un); 
    }
    printf("%d\n",qMax);
}
double m11(double deltaT, double U1)
{
    double alfa = beta * N - lambda;
    return 1 - deltaT * a11 * dfun(U1);
}
double m12(double deltaT, double U2)
{
    double alfa = beta * N - lambda;
    return (-1) * deltaT * a12 * dfun(U2);
}
double m21(double deltaT, double U1)
{
    double alfa = beta * N - lambda;
    return (-1) * deltaT * a21 * dfun(U1);
}
double m22(double deltaT, double U2)
{
    double alfa = beta * N - lambda;
    return 1 - deltaT * a22 * dfun(U2);
}
double F1(double U1, double U2, double un, double deltaT)
{
    return U1 - un - deltaT * (a11 * fun(U1) + a12 * fun(U2));
}
double F2(double U1, double U2, double un, double deltaT)
{
    return U2 - un - deltaT * (a21 * fun(U1) + a22 * fun(U2));
}
double deltaU1(double U1, double U2, double un, double deltaT)
{
    return (F2(U1,U2,un,deltaT) * m12(deltaT,U2) - F1(U1,U2,un,deltaT) * m22(deltaT,U2))
    /
    (m11(deltaT,U1) * m22(deltaT, U2) - m12(deltaT,U2) * m21(deltaT,U1));
}
double deltaU2(double U1, double U2, double un, double deltaT)
{
    return (F1(U1,U2,un,deltaT) * m21(deltaT,U1) - F2(U1,U2,un,deltaT) * m11(deltaT,U1))
    /
    (m11(deltaT,U1) * m22(deltaT, U2) - m12(deltaT,U2) * m21(deltaT,U1));
}
void RK2(FILE * fp)
{
    double tmax = 100;
    double deltaT = 0.1;
    
    int q  = 0;
    int qMax = 0;  

    double b1 = 0.5;
    double b2 = b1;

    double un = 1.;
    double U1 = 1.;
    double U2 = 1.;
    double U1q = 0.;
    double U2q = 0.;
    for(double t = deltaT; t <= tmax; t += deltaT)
    {
        U1q = 0.;
        U2q = 0.;
        while((fabs(U1 - U1q) > TOL || fabs(U2 - U2q) > TOL) && q < 20  )
        {
            U1q = U1;
            U2q = U2;
            U1 = U1q + deltaU1(U1,U2,un,deltaT);
            U2 = U2q + deltaU2(U1,U2,un,deltaT);
            q++;
        }
        un = un + deltaT*( b1 * fun(U1) + b2 * fun(U2));
        fprintf(fp,"%g %g %g %d\n",t,un, N - un,q); 
        if(q > qMax)  
            qMax = q;
        q=0;
    }
    printf("%d\n",qMax);
}
int main(void) {
    FILE * fp = fopen("out.txt","w");
	MetodaPicarda(fp);
    fclose(fp);
    fp = fopen("out2.txt","w");
    MetodaNewtona(fp);
    fclose(fp);
    fp = fopen("out3.txt","w");
    RK2(fp);
    fclose(fp);
	return 0;
}
