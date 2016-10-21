#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define tMAX 20.
#define DELTA 0.0005

const int n = 8;
const int N = n*n*n;           //total number of particles
double MASS_ = 1.;             //mass
double initV2 = 0.5 , initV = sqrt(initV2);      //скорость теплового движени€
const double rmin = pow(10, (1. / 3));          //посто€нна€ решетки
const double l_ = n * rmin;                //size of cell

const double dr = 0.01, V_tot = l_*l_*l_;         //шаг дл€ опред бин коррел€ции
const int N_r = (int)(l_ / 2. / dr);            //размер массива дл€ бинарной коррел€ции
double lambda, deltE, E_av, E_= 3./ 2. * 1.3;
double E_os = 0.8 / sqrt((double)N);              //допустимые колебани€ энергии (дл€ выключени€ термостата)
double dt = (double)DELTA;
int Nt = (int)(tMAX/DELTA);
int rel = (int)(10./dt);             //врем€ ожидани€ термостата (если колебани€ энергии < E_os в течение 10, то термостат выключаетс€)
double Utot, Ekin;
double r[N][3], v[N][3], F[N][3], m[N], r_n[N][3], L[3], L_2[3];	
double g[N_r];
FILE *f, *f_b;

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
	fprintf(f,"%lf %lf\n", Utot, Ekin);
}

inline void save_bin_cor()
{
	int i;
	double r, V_k = 0.01*4./3.*3.14159265*(l_/2.)*(l_/2.)*(l_/2.);
	for (i=0;i<N_r;i++)
	{
		r = (i + 0.5)*dr;
		fprintf(f_b,"%lf %lf\n", r, (g[i] / (4.*3.14159265*r*r*dr*Nt*V_k)) );
	}
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
	int i, j, k, count = 0;
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

inline char termo_B()
{
	E_av = Ekin * 1. / N;
	deltE = fabs(E_av - E_) * 1. / E_;
	int i,j;
	lambda = sqrt(1. + (0.0005/1.) * (E_ * 1. / E_av - 1));
	for (i = 1; i < N; i++)
		for (j=0;j<3;j++)
		{
			v[i][j] *= lambda;
		}
	return (deltE < E_os);		
}

inline void clearBin()
{
	int i;
	for (i=0;i<N_r;i++)
	{
		g[i]=0;
	}
}

inline void bin_cor()
{
	int i,j,k,t;
	double r_v[3], r2;
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
			if ((t = (int)(floor(sqrt(r2)/dr))) < N_r)
				g[t]+=1;
		}
}
  
int main(void)
{
	srand(time(NULL));	
	int i,j,k,n,t,counter;
	for (i=0;i<3;i++)
	{
		L[i] = l_;
		L_2[i] = l_ / 2.;
	}
	for (i=0;i<N;i++)
	{
		m[i]=MASS_;
	}
	defPos_Cryst();
	defVel();	
	f=fopen("statistic.txt", "w");
	f_b=fopen("bin_cor.txt", "w");
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
	saveEnergy();          //начальна€ энерги€
	while (1)          //цикл с термостатом и условием выключени€
	{		
		clearF();
		n_image();
		calcF();
		eqMot();
		saveEnergy();
		if (termo_B())
		{
			counter+=1;
			if (counter == rel)
				break;
		}
		else
			counter = 0;
	}
	printf("general loop\n");
	clearBin();
	for (n=0; n<Nt; n++)          //основной цикл
	{	
		bin_cor();	
		clearF();
		n_image();
		calcF();
		eqMot();
		saveEnergy();
	}	
	save_bin_cor();
	t = clock() - t;
	printf("%f\n", (float)(t * 1. / CLOCKS_PER_SEC));	
	fclose(f);
	fclose(f_b);
	system("pause");
	return 0;
}
