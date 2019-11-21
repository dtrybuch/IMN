#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double firstFunction(double lambda,double y)
{
    return lambda*y;
}
double firstRealFunction(double lambda,double t)
{
    return exp(lambda*t);
}
double EulerFunction(double deltaT,double y,double lambda)
{
    return y + deltaT * firstFunction(lambda,y);
}

void FirstExercise(double deltaT)
{
    double lambda = -1;
    double y = 1;
    FILE *fp=fopen("out.dat", "w");
    for(double x = 0.; x <= 5.; x = x + deltaT)
    {
        y =  EulerFunction(deltaT,y,lambda);
        fprintf(fp,"%g ",x + deltaT);
        fprintf(fp,"%g ",y);
        fprintf(fp,"%g\n",y - firstRealFunction(lambda,x+deltaT));            
    }
    fprintf(fp,"\n");
    fclose (fp);
}
void RK2(double deltaT)
{
    double lambda = -1;
    double y = 1;
    double k1,k2;
    FILE *fp=fopen("out2.dat", "w");
    for(double x = 0.; x <= 5.; x = x + deltaT)
    {
        k1 = firstFunction(lambda,y);
        k2 = firstFunction(lambda,y + deltaT*k1);
        y = y + (deltaT/2) * (k1+k2);
        fprintf(fp,"%g ",x + deltaT);
        fprintf(fp,"%g ",y);
        fprintf(fp,"%g\n",y - firstRealFunction(lambda,x+deltaT));            
    }
    fprintf(fp,"\n");
    fclose (fp);
}
void RK4(double deltaT)
{
    double lambda = -1;
    double y = 1;
    double k1,k2,k3,k4;
    FILE *fp=fopen("out3.dat" , "w");
    for(double x = 0.; x <= 5.; x = x + deltaT)
    {
        k1 = firstFunction(lambda,y);
        k2 = firstFunction(lambda,y + (deltaT/2)*k1);
        k3 = firstFunction(lambda,y + (deltaT/2)*k2);
        k4 = firstFunction(lambda,y + (deltaT*k3));
        y = y + (deltaT/6.)*(k1 + 2 * k2 + 2 * k3 + k4);
        fprintf(fp,"%g ",x + deltaT);
        fprintf(fp,"%g ",y);
        fprintf(fp,"%g\n",y - firstRealFunction(lambda,x+deltaT));            
    }
    fprintf(fp,"\n");
    fclose (fp);
}
double VFunction(double omegaV,double t)
{
    return 10*sin(omegaV*t);
}
double GFunction(double t,double L,double R, double I, double Q, double omegaV,double C )
{
    return VFunction(omegaV,t)/L - (R/L)*I - (1./(L*C))*Q;
}
double IFunction(double I, double deltaT,double k)
{
    return I + deltaT*k;
}
void RRZ2(double mnoznik)
{
    double k1q,k1i,k2q,k2i,k3q,k3i,k4q,k4i;
    double deltaT = 0.0001;
    double R = 100;;
    double L = 0.1;
    double C = 0.001;
    double omega = 1./(sqrt(L*C));
    double omegaV = mnoznik*omega;
    double T0 = (2. * M_PI)/omega; 
    double Q = 0.;
    double I = 0.;
    FILE *fp=fopen("out4.dat" , "w");
    for(double x = 0.; x <= 4*T0; x = x + deltaT)
    {
        k1q = I;
        k1i = GFunction(x + deltaT,L,R,I,Q,omegaV,C);
        k2q = IFunction(I,deltaT/2.,k1i);
        k2i = GFunction(x + (1/2.)*deltaT,L,R,I + (deltaT/2.)*k1i,Q + (deltaT/2.)*k1q,omegaV,C);
        k3q = IFunction(I,deltaT/2.,k2i);
        k3i = GFunction(x + (1/2.)*deltaT,L,R,I + (deltaT/2.)*k2i,Q + (deltaT/2.)*k2q,omegaV,C);
        k4q = IFunction(I,deltaT,k3i);
        k4i = GFunction(x + deltaT,L,R,I + deltaT*k3i,Q + deltaT*k3q,omegaV,C);
        Q = Q + (deltaT/6.)*(k1q + 2*k2q + 2*k3q + 2*k4q);
        I = I + (deltaT/6.)*(k1i + 2*k2i + 2*k3i + 2*k4i);
        fprintf(fp,"%g ",x + deltaT);
        fprintf(fp,"%g ",I);
        fprintf(fp,"%g\n",Q);           
    }
    fprintf(fp,"\n");
    fclose (fp);
}
int main(void) {
	//metoda Eulera
    FirstExercise(0.01);
    FirstExercise(0.1);
    FirstExercise(1.0);
    //metoda trapezow
    RK2(0.01);
    RK2(0.1);
    RK2(1.);
    //metoda jawana RK4
    RK4(0.01);
    RK4(0.1);
    RK4(1.);
    //RRZ 2 rzedu
    RRZ2(0.5);
    RRZ2(0.8);
    RRZ2(1.0);
    RRZ2(1.2);

	return 0;
}
