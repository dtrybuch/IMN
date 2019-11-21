#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tmax 40
#define gamma 1e-10

double max(double x, double y){
    return x > y ? x : y;
}
double f(double v)
{
	return v;
}
double g(double x, double v, double alfa)
{
	return alfa * (1 - x * x) * v - x;
}
double F(double xnext, double xn, double deltaT, double vn, double vnext, double alfa)
{
	return xnext - xn - (deltaT/2.) * (f(vn) + f(vnext));
}
double G(double xnext, double xn, double deltaT, double vn, double vnext,double alfa)
{
	return vnext - vn - (deltaT/2.) * (g(xn,vn,alfa) + g(xnext,vnext,alfa));
}

double a11(){
    return 1.;
}

double a12(double deltaT){
    return (-1)*(deltaT/2.);
}

double a21(double deltaT, double xnext, double vnext, double alfa){
    return (-1)*( deltaT / 2. )*((-2) * alfa * xnext * vnext - 1);
}

double a22(double deltaT, double xnext, double alfa){
    return (1 - (deltaT/2.)*alfa*(1-xnext * xnext));
}
void MetodaTrapezow(double *xn, double *vn, double deltaT, double alfa)
{
	double xk = *xn;
	double vk = *vn;
	double xnext = *xn;
	double vnext = *vn;

	double deltaX = 1000;
	double deltaV = 10000;
	while(fabs(deltaX) > gamma || fabs(deltaV) > gamma)
	{
		deltaX = 0.;
		deltaV = 0.;

		deltaX = (-F(xnext,*xn,deltaT,*vn,vnext,alfa)*a22(deltaT,xnext,alfa) - (-G(xnext,*xn,deltaT,*vn,vnext,alfa)*a12(deltaT)))
		/
		(a11()*a22(deltaT,xnext,alfa) - a12(deltaT)*a21(deltaT,xnext,vnext,alfa));

		deltaV = (a11()*(-G(xnext,*xn,deltaT,*vn,vnext,alfa) - a21(deltaT,xnext,vnext,alfa) * (-F(xnext,*xn,deltaT,*vn,vnext,alfa))))
		/
		(a11()*a22(deltaT,xnext,alfa) - a12(deltaT)*a21(deltaT,xnext,vnext,alfa));
		xnext = xnext + deltaX;
		vnext = vnext + deltaV;
	}
	*xn = xnext;
	*vn = vnext;
}

void RK2(double *xn, double *vn, double deltaT, double alfa)
{
	double k1x = f(*vn);
	double k1v = g(*xn,*vn,alfa);

	double k2x = f(*vn + deltaT * k1v);
	double k2v = g(*xn + deltaT * k1x, *vn + deltaT * k1v, alfa);
	*xn = *xn + (deltaT/2.) * (k1x + k2x);
	*vn = *vn + (deltaT/2.) * (k1v + k2v);
}

void KontrolaKroku(FILE * fp,double TOL, void(*schemat)(double *xn, double *vn, double deltaT, double alfa),double alfa)
{
	double Ex = 0;
	double Ev = 0;
	double x0 = 0.01;
	double v0 = 0;
	double xn = x0;
    double vn = v0;
	double S = 0.75;
	double p = 2.;
	double deltaT0 = 1.;
	double t = 0;
	double deltaT = deltaT0;
	double xtmp = 0;
	double vtmp = 0;
	double xtmp2 = 0;
	double vtmp2 = 0;
	fprintf(fp, "%f %f %f %f \n", t, deltaT, xn, vn);
	do
	{
		xtmp = xn;
		vtmp = vn;
		xtmp2 = xn;
		vtmp2 = vn;
		schemat(&xtmp,&vtmp,deltaT,alfa);
		schemat(&xtmp,&vtmp,deltaT,alfa);
		schemat(&xtmp2, &vtmp2, 2 * deltaT, alfa);
		Ex = (xtmp - xtmp2) / (pow(2,p) - 1);
		Ev = (vtmp - vtmp2) / (pow(2,p) - 1);
		if(max(fabs(Ex),fabs(Ev)) < TOL ){
            t = t + 2 * deltaT;
            xn = xtmp;
            vn = vtmp;
            fprintf(fp, "%f %f %f %f \n", t, deltaT, xn, vn);
		}
		deltaT = (pow((S * TOL) / (max(fabs(Ex),fabs(Ev))), (1./(p+1)))) * deltaT;
	} while (t < tmax );
}
int main(void) {
	FILE * fp = fopen("RK21.txt","w"); 
	
	double TOL = pow(10,-2);
	double alfa = 5;
// zadanie 1
    KontrolaKroku(fp, TOL, RK2,alfa);
	fclose(fp);

	fp = fopen("RK22.txt","w");
    TOL = pow(10,-5);
    KontrolaKroku(fp, TOL, RK2, alfa);
	fclose(fp);

	fp = fopen("trapezy1.txt", "w");
    TOL = pow(10,-2);
    KontrolaKroku(fp, TOL, MetodaTrapezow,alfa);
	fclose(fp);

	fp = fopen("trapezy2.txt", "w");
    TOL = pow(10,-5);
    KontrolaKroku(fp, TOL, MetodaTrapezow,alfa);
	fclose(fp);
	return 0;
}
