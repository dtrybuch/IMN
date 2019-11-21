#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-8
#define epsilon 1
#define delta 0.1
#define nx 150
#define ny 100
#define V1 10.
#define V2 0.
#define xmax (delta * nx) 
#define ymax (delta * ny)

double ro[nx+1][ny+1];
void gestosc()
{
	double deltax = 0.1 * xmax;
	double deltay = 0.1 * ymax;
	for(int i = 0;i <= nx;i++)
	{
		for (int j = 0;j <=ny;j++)
		{
			ro[i][j] = exp((-1)* pow(delta*i - 0.35 * xmax, 2) / pow(deltax,2) - pow(delta*j - 0.5 * ymax,2) / pow(deltay,2))
					- exp((-1) * pow(delta*i - 0.65 * xmax, 2) / pow(deltax,2) - pow(delta*j - 0.5 * ymax,2) / pow(deltay,2));
		}
	}
}

double stop(double Vn[nx + 1][ny+1])
{
	
	double s = 0.;
	for(int i = 0; i < nx; i++)
	{
		for(int j = 0; j < ny; j++)
		{
			s += pow(delta,2) * (1./2. * (pow( (Vn[i+1][j] - Vn[i][j]) / delta,2))
			  	+ 1./2. * (pow( (Vn[i][j+1] - Vn[i][j]) / delta,2)) 
			  	- ro[i][j] * Vn[i][j] );
		}
	}
	return s;
}
double blad(FILE * fp, double V[nx+1][ny+1])
{

	double result = 0;
	for(int i = 1; i< nx; i++)
	{
		for(int j = 1;j < ny; j++)
		{
			result = (V[i+1][j] - 2*V[i][j] + V[i-1][j])/pow(delta,2) + (V[i][j+1] - 2*V[i][j] + V[i][j-1])/pow(delta,2) + ro[i][j]/epsilon;
			fprintf(fp,"%g %g %g\n",i*delta,j*delta,result);
		}
		fprintf(fp,"\n");
	}
}
void globalna(FILE * fp,FILE * fp2,FILE * fp3,double omegaG)
{
	double Vn[nx+1][ny+1];
	double Vs[nx+1][ny+1];
	int it = 0;
	double s = 1.;
	double snext = 0.;
	for(int i = 0;i <= nx; i++)
	{
		Vn[i][0] = V1;
		Vn[i][ny] = V2;
		Vs[i][0] = V1;
		Vs[i][ny] = V2;
	}

	for(int i = 0; i <= nx; i++)
	{
		for(int j = 1; j < ny; j++)
		{
			Vn[i][j] = 0.;
			Vs[i][j] = 0.;
		}
	}
	while(fabs((snext - s) / s) > TOL)
	{
		s = snext;
		for(int i = 1;i < nx; i++)
		{ 
			for(int j = 1; j < ny; j++)
			{
				Vn[i][j] = (1./4.) * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + pow(delta,2) / epsilon * ro[i][j]) ; 
			}
		}
		for(int j = 1; j < ny; j++)
		{
			Vn[0][j] = Vn[1][j];
			Vn[nx][j] = Vn[nx-1][j];
		}

		for(int i = 0; i <= nx; i++)
		{
			for(int j = 0; j <= ny; j++)
			{
				Vs[i][j] = (1 - omegaG) * Vs[i][j] + omegaG * Vn[i][j];
			}
		}
		snext = stop(Vs);
		it++;
		fprintf(fp,"%d %g \n",it, snext);
		//printf("%d %g \n",it, snext);
		// printf("     %g \n", fabs((snext - s) / s));
	}
	for(int i = 0; i<=nx; i++)
	{
		for(int j = 0; j<= ny; j++)
		{
			fprintf(fp2,"%g %g %g \n",i*delta,j*delta,Vs[i][j]);
		}
		fprintf(fp2,"\n");
	}
	blad(fp3,Vs);
	printf("%d\n",it);
}
void lokalna(FILE * fp,FILE * fp2, FILE * fp3,double omegaL)
{
	double Vn[nx+1][ny+1];
	int it = 0;
	double s = 1.;
	double snext = 0.;
	for(int i = 0;i <= nx; i++)
	{
		Vn[i][0] = V1;
		Vn[i][ny] = V2;
	}

	for(int i = 0; i <= nx; i++)
	{
		for(int j = 1; j < ny; j++)
		{
			Vn[i][j] = 0.;
		}
	}
	while(fabs((snext - s) / s) > TOL)
	{
		s = snext;
		for(int i = 1;i < nx; i++)
		{ 
			for(int j = 1; j < ny; j++)
			{
				Vn[i][j] = (1 - omegaL) * Vn[i][j] + (omegaL/4.) * (Vn[i+1][j] + Vn[i-1][j] + Vn[i][j+1] + Vn[i][j-1] + (pow(delta,2) / epsilon) * ro[i][j]);
			}
		}
		for(int j = 1; j < ny; j++)
		{
			Vn[0][j] = Vn[1][j];
			Vn[nx][j] = Vn[nx-1][j];
		}
		snext = stop(Vn);
		it++;
		fprintf(fp,"%d %g \n",it, snext);
	}
	printf("%d\n",it);
}
int main(void) {
	gestosc();
	FILE * fp, *fp2, *fp3;
	double omegaG;
	double omegaL;
	fp = fopen("out.txt","w");
	fp2 = fopen("outV.txt","w");
	fp3 = fopen("outBlad.txt","w");
	omegaG = 0.6;
	globalna(fp,fp2,fp3,omegaG);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);

	fp = fopen("out2.txt","w");
	fp2 = fopen("outV2.txt","w");
	fp3 = fopen("outBlad2.txt","w");
	omegaG = 1.0;
	globalna(fp,fp2,fp3,omegaG);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);

	fp = fopen("out3.txt","w");
	fp2 = fopen("outV3.txt","w");
	fp3 = fopen("outBlad3.txt","w");
	omegaL = 1.0;
	lokalna(fp,fp2,fp3,omegaL);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);

	fp = fopen("out4.txt","w");
	fp2 = fopen("outV4.txt","w");
	fp3 = fopen("outBlad4.txt","w");
	omegaL = 1.4;
	lokalna(fp,fp2,fp3,omegaL);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);

	fp = fopen("out5.txt","w");
	fp2 = fopen("outV5.txt","w");
	fp3 = fopen("outBlad5.txt","w");
	omegaL = 1.8;
	lokalna(fp,fp2,fp3,omegaL);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);

	fp = fopen("out6.txt","w");
	fp2 = fopen("outV6.txt","w");
	fp3 = fopen("outBlad6.txt","w");
	omegaL = 1.9;
	lokalna(fp,fp2,fp3,omegaL);
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	return 0;
}
