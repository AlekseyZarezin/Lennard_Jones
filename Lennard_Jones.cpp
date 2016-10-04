#include<stdio.h>
#include<stdlib.h>
//#include<string.h>
#include<math.h>
#include<time.h>

#define tMAX 10.
#define DELTA 0.0005
#define N 4
#define EPS 1.
#define SIGMA 1.
#define MASS_ 1.

double dt = (double)DELTA;
int Nt = (int)(tMAX/DELTA);
double Utot, Ekin;
double r[N][3], v[N][3], F[N][3], m[N];	
double sigm6 = (double)(SIGMA*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA);
FILE *f;	

inline double U(double r2)                                        //вычисление потенциала
{
	if (r2<=(9*SIGMA*SIGMA))
	{
		double r6 = r2*r2*r2;
		return 4*(double)EPS*((sigm6*sigm6/(r6*r6)) - (sigm6/r6));
	}
	else
	{
		return 0.;
	}
}

double F_r(double r2)
{
	if (r2<=(9*SIGMA*SIGMA))
	{
		double r6 = r2*r2*r2;		
		return 4*(double)EPS*( 12.*(sigm6*sigm6/(r6*r6*r2)) - 6.*(sigm6/(r6*r2)) );
	}
	else
	{
		return 0.;
	}
}

void clearF()
{
	int i,j;
	for (i=0;i<N;i++)
		for (j=0;j<3;j++)
		{
			F[i][j]=0;
		}
}

void calcF()
{
	int i,j,k;
	Utot = 0;
	double r_v[3], r2, f_r;
	for (i=0;i<N-1;i++)
		for (j=i+1;j<N;j++)
		{
			r2=0;
			for (k=0;k<3;k++)
			{
				r_v[k] = r[i][k] - r[j][k];
				r2 += r_v[k] * r_v[k];
			}
			f_r = F_r(r2);
			Utot += U(r2);
			for (k=0;k<3;k++)
			{
				F[i][k] += f_r*r_v[k];
				F[j][k] -= f_r*r_v[k];
			}
		}
}

void calcEkin()
{
	double v2;
	int i,k;
	Ekin = 0;
	for (i=0;i<N;i++)
	{
		v2 = 0;
		for (k=0;k<3;k++)
			v2 += v[i][k] * v[i][k];
		Ekin += m[i] * v2 / 2.;
	}
}

void eqMot()
{
	int i,k;
	double v2;
	Ekin = 0;
	for (i=0;i<N;i++)
	{
		v2 = 0;
		for (k=0;k<3;k++)
		{
			v2 += (v[i][k] + 0.5 * F[i][k] * dt / m[i]) * (v[i][k] + 0.5 * F[i][k] * dt / m[i]);
			v[i][k]+=F[i][k] * dt / m[i];
			r[i][k]+=v[i][k] * dt;
		}
		Ekin += v2 * m[i] /2.;
	}
}

void saveEnergy()
{
	f=fopen("statistic.txt", "a");
	fprintf(f,"%lf %lf\n", Utot, Ekin);
	fclose(f);
}

int main(void)
{	
	int i,j,k,n;
	for (i=0;i<N;i++)
	{
		m[i]=(double)MASS_;
	}
	for (i=0;i<N;i++)
	{
		scanf("%lf %lf %lf", &(r[i][0]), &(r[i][1]), &(r[i][3]));
	}
	for (i=0;i<N;i++)
	{
		scanf("%lf %lf %lf", &(v[i][0]), &(v[i][1]), &(v[i][3]));
	}
	
	//printf("%lf\n", U(0.9));
	
	calcEkin();                   
	clearF();
	calcF();
	for (i=0;i<N;i++)                   //значение скорости V(1/2)
		for (k=0;k<3;k++)
		{
			v[i][k] += 0.5 * F[i][k] * dt / m[i];
			r[i][k] += v[i][k] * dt;
		}
	saveEnergy();          //начальная энергия
		
	for (n=0; n<Nt; n++)          //основной цикл
	{		
		clearF();
		calcF();
		//calcEkin();
		//save v[][],r[][],F[][] in file
		eqMot();
		saveEnergy();
	}
	
	system("pause");
	return 0;
}
