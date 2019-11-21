#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-8
#define epsilon 1
#define delta 0.2
#define nx 128
#define ny 128
#define xmax (delta * nx) 
#define ymax (delta * ny)
double VB13(double y)
{
	return sin(M_PI * y / ymax);
}
double VB2(double x)
{
	return  - sin(2 * M_PI * x / xmax);
}
double VB4(double x)
{
	return  sin(2 * M_PI * x / xmax);
}
double stop(double V[nx + 1][ny + 1], int k)
{
	double s = 0;
	for( int i = 0; i <= nx - k; i+=k)
	{
		for(int j = 0; j <= ny - k; j+=k)
		{
			s +=  (pow(k * delta,2) / 2.) 
			* 
			(
			pow((V[i+k][j] - V[i][j]) / (2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2 * k * delta),2)   
			+
			pow((V[i][j+k] - V[i][j]) / ( 2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta),2)
			);
		}
	}
	return s;
}
void Relaksacja(FILE *fpV, FILE * fpM, int k)
{
	double V[nx + 1][ny + 1];
	double S = 0;
	double Snext = 1;
	int iter = 0;
	for(int i = 0; i <= nx; i++)
	{
		V[i][ny] = VB2(delta * i);
		V[i][0] = VB4(delta * i);
	}
	for(int j = 0; j <= ny; j++)
	{
		V[0][j] = VB13(delta * j);
		V[nx][j] = VB13(delta * j);
	}
	for(int i = 1; i < nx; i++)
	{
		for(int j = 1; j < ny; j++)
		{
			V[i][j] = 0;
		}
	}
	while(k >= 1)
	{
		Snext = 5316131;
		while(fabs((Snext - S) / S) > TOL)
		{
			Snext = S;
			for (int i = k; i <= nx - k; i+=k )
			{
				for (int j = k; j <= ny - k; j+=k )
				{
					V[i][j] = 0.25 * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
				}
			}
			S = stop(V,k);
			iter++;
			fprintf(fpV,"%d %d %f \n",k,iter,S);
		}

		fprintf(fpV,"\n\n");

		for(int i = 0; i <=nx; i+=k)
		{
			for(int j = 0; j<=ny; j+=k)
			{
				fprintf(fpM,"%f %f %f \n", delta*i, delta*j, V[i][j]);
			}
			fprintf(fpM,"\n");
		}
		fprintf(fpM,"\n");

		for(int i = 0; i<= nx - k; i+=k)
		{
			for(int j = 0; j <= ny - k; j+=k)
			{
				V[i + k/2][j + k/2] = 0.25 * (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
				if(i != 0)
					V[i][j + k/2] = 0.5 * (V[i][j] + V[i][j+k]);
				if(i != nx-k)
					V[i + k][j + k/2] = 0.5 * (V[i+k][j] + V[i+k][j+k]);
				if(j != 0)
					V[i + k/2][j] = 0.5 * (V[i][j] + V[i+k][j]);
				if(j != ny-k)
					V[i+k/2][j+k] = 0.5 *(V[i][j+k] + V[i+k][j+k]);
			}
		}
		k = k/2;
	}
}
int main(void) {
	int k = 16;
	FILE * fpV = fopen("V.dat", "w");
	FILE * fpM = fopen("map.dat","w");
	Relaksacja(fpV, fpM, k);
	fclose(fpM);
	fclose(fpV);
	return 0;
}
