#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define tMAX 10.
#define DELTA 0.0005

const int N = 100;
double MASS_ = 1.;             //mass
double l_ = 10.;                //size of cell
double initV2 = 3.* 1.3 / MASS_, initV = sqrt(initV2);
double rmin = pow(10, (1. / 3));          //постоянная решетки

double dt = (double)DELTA;
int Nt = (int)(tMAX/DELTA);
double Utot, Ekin;
double r[N][3], v[N][3], F[N][3], m[N], r_n[N][3], L[3], L_2[3];	
FILE *f;

inline double U(double r2)                                        //вычисление потенциала
{
	if (r2<=9.)
	{
		double r6 = r2*r2*r2;
		return 4.*((1./(r6*r6)) - (1./r6)) + 0.00547944;
	}
	else
	{
		return 0.;
	}
}

inline double F_r(double r2)
{
	if (r2<=9.)
	{
		double r6 = r2*r2*r2;		
		return 4.*( (12./(r6*r6*r2)) - (6./(r6*r2)) );
	}
	else
	{
		return 0.;
	}
}

inline void n_image()
{
	int i,j;
	for (i=0;i<N;i++)
		for (j=0;j<3;j++)
		{
			if (r[i][j]>0)
				r_n[i][j] = fmod(r[i][j] + L_2[j], L[j]) - L_2[j];
			else
				r_n[i][j] = fmod(r[i][j] - L_2[j], L[j]) + L_2[j];
		}
}

inline void clearF()
{
	int i,j;
	for (i=0;i<N;i++)
		for (j=0;j<3;j++)
		{
			F[i][j]=0;
		}
}

inline void calcF()
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
				r_v[k] = r_n[i][k] - r_n[j][k];
				if (r_v[k] > L_2[k])
					r_v[k] -= L[k];
				else
					if (r_v[k] < (-1.)*L_2[k])
						r_v[k] += L[k];
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

inline void calcEkin()
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

inline void eqMot()
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

inline void saveEnergy()
{
	fprintf(f,"%lf %lf\n", Utot, Ekin, r[1][1]);
}

inline void defVel()
{
	int i;	
	for (i=0;i<N;i++)
	{
		v[i][0] = (rand() * 2. / RAND_MAX - 1.) * initV;
		v[i][1] = (rand() * 2. / RAND_MAX - 1.) * sqrt(initV2 - v[i][0] * v[i][0]);
		v[i][2] = sqrt( fabs(initV2 - v[i][0] * v[i][0] - v[i][1] * v[i][1]) );
	}
}

inline void defPos_Cryst()
{
	int i, j, k, count = 0, n = (int)pow((double)(N), (1./3))+2;
	printf("%d\n", n);
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			for (k=0;k<n;k++)
			{
				r[count][0] = (double)(i*rmin);
				r[count][1] = (double)(j*rmin);
				r[count][2] = (double)(k*rmin);
				count++;
				if (count==N)
					return;
			}
}
  
int main(void)
{
	srand(time(NULL));	
	int i,j,k,n,t;
	for (i=0;i<3;i++)
	{
		L[i] = l_;
		L_2[i] = l_ / 2.;
	}
	for (i=0;i<N;i++)
	{
		m[i]=MASS_;
	}
	//for (i=0;i<N;i++)
	//{
	//	scanf("%lf %lf %lf", &(r[i][0]), &(r[i][1]), &(r[i][3]));
	//}
	//for (i=0;i<N;i++)
	//{
	//	scanf("%lf %lf %lf", &(v[i][0]), &(v[i][1]), &(v[i][3]));
	//}
	defPos_Cryst();
	defVel();
	
	f=fopen("statistic.txt", "w");

	t = clock();	
	calcEkin();                   
	clearF();
	n_image();
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
		n_image();
		calcF();
		//calcEkin();
		//save v[][],r[][],F[][] in file
		eqMot();
		saveEnergy();
	}
	t = clock() - t;
	printf("%f\n", (float)(t * 1. / CLOCKS_PER_SEC));
	
	fclose(f);
	system("pause");
	return 0;
}
