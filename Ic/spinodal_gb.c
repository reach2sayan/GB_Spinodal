/************************************************************h
To compile
gcc "binary.c" -LLIBDIR=/usr/local/lib -lgsl -lgslcblas -lfftw3 -lm -o "ternary.out"

Author: Sayan Samanta

************************************************************/

//variables to take input X, Ka, Kb, cA, cB, Maa, Mbb, Mab	

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>

int main(void) 
{
	FILE *gp, *fpi, *fpo, *fpon1, *fpon2;
	
	gp = fopen("plotAnimation","w");
	fpi = fopen("inputFile.dat","r");
	
	fprintf(gp, "set cbrange[0:1]\n");
	fprintf(gp, "set palette defined (0 'black', 1 'white')\n");

	int i1,i2;
	int half_nx, half_ny;
	double kx, ky, delta_kx, delta_ky;
	double kx2, ky2, k2, k4;
	char name[50];
	int INDEX, TPRINT = 100;

	double A = 1.0, kc = 1.0, M = 1.0, e = 2.0, kn = 0.5, L = 1.0;
	double delta_x;
	double delta_y;
	double delta_t;
	double icA = 0.5, icB = 0.5;
	delta_x = 1.0; delta_y = 1.0; delta_t = 0.1;
	int n_x, n_y;

	n_x = 512;
	n_y = 256;

	fftw_complex *c, *n1, *n2, *h1, *h2, *g, *ctemp, *n1temp, *n2temp;
	fftw_plan planFc, planBc, planFn1, planBn1, planFn2, planBn2, planFg, planBg, planFh1, planBh1, planFh2, planBh2, planFctemp, planBctemp, planFn1temp,planBn1temp, planFn2temp, planBn2temp;

	c = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	n1 = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	n2 = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	h1 = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	h2 = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	g = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	ctemp = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	n1temp = fftw_malloc(n_x*n_y*sizeof(fftw_complex));
	n2temp = fftw_malloc(n_x*n_y*sizeof(fftw_complex));

	planFc = fftw_plan_dft_2d(n_x,n_y,c,c,FFTW_FORWARD,FFTW_ESTIMATE);
	planFn1 = fftw_plan_dft_2d(n_x,n_y,n1,n1,FFTW_FORWARD,FFTW_ESTIMATE);
	planFn2 = fftw_plan_dft_2d(n_x,n_y,n2,n2,FFTW_FORWARD,FFTW_ESTIMATE);
	planFg = fftw_plan_dft_2d(n_x,n_y,g,g,FFTW_FORWARD,FFTW_ESTIMATE);
	planFh1 = fftw_plan_dft_2d(n_x,n_y,h1,h1,FFTW_FORWARD,FFTW_ESTIMATE);
	planFh2 = fftw_plan_dft_2d(n_x,n_y,h2,h2,FFTW_FORWARD,FFTW_ESTIMATE);
	planFctemp = fftw_plan_dft_2d(n_x,n_y,ctemp,ctemp,FFTW_FORWARD,FFTW_ESTIMATE);
	planFn1temp = fftw_plan_dft_2d(n_x,n_y,n1temp,n1temp,FFTW_FORWARD,FFTW_ESTIMATE);
	planFn2temp = fftw_plan_dft_2d(n_x,n_y,n2temp,n2temp,FFTW_FORWARD,FFTW_ESTIMATE);
	
	planBc = fftw_plan_dft_2d(n_x,n_y,c,c,FFTW_BACKWARD,FFTW_ESTIMATE);
	planBn1 = fftw_plan_dft_2d(n_x,n_y,n1,n1,FFTW_BACKWARD,FFTW_ESTIMATE);
	planBn2 = fftw_plan_dft_2d(n_x,n_y,n2,n2,FFTW_BACKWARD,FFTW_ESTIMATE);
	planBg = fftw_plan_dft_2d(n_x,n_y,g,g,FFTW_BACKWARD,FFTW_ESTIMATE);
	planBh1 = fftw_plan_dft_2d(n_x,n_y,h1,h1,FFTW_BACKWARD,FFTW_ESTIMATE);
	planBh2 = fftw_plan_dft_2d(n_x,n_y,h2,h2,FFTW_BACKWARD,FFTW_ESTIMATE);
	planBctemp = fftw_plan_dft_2d(n_x,n_y,ctemp,ctemp,FFTW_FORWARD,FFTW_ESTIMATE);
	planBn1temp = fftw_plan_dft_2d(n_x,n_y,n1temp,n1temp,FFTW_FORWARD,FFTW_ESTIMATE);
	planBn2temp = fftw_plan_dft_2d(n_x,n_y,n2temp,n2temp,FFTW_FORWARD,FFTW_ESTIMATE);
	
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for(i1=0; i1 < n_x; ++i1)
	{
		for(i2=0; i2 < n_y; ++i2)
		{
			double u = gsl_rng_uniform(r);
			__real__ c[i2+n_y*i1] = icA + (icA - u)*5e-3;
			__imag__ c[i2+n_y*i1] = 0.0;
					
			__real__ ctemp[i2+n_y*i1] = 0.0;
			__imag__ ctemp[i2+n_y*i1] = 0.0;
			

			if(i1 < n_x/2)
			  {
			    __real__ n1[i2 + n_y*i1] = 0.0;
			    __real__ n2[i2 + n_y*i1] = 1.0;
			    __imag__ n1[i2 + n_y*i1] = 0.0;
			    __imag__ n2[i2 + n_y*i1] = 0.0;
			  }
			else
			  {
			    __real__ n1[i2 + n_y*i1] = 1.0;
			    __real__ n2[i2 + n_y*i1] = 0.0;
			    __imag__ n1[i2 + n_y*i1] = 0.0;
			    __imag__ n2[i2 + n_y*i1] = 0.0;
			  }

			//__real__ m[i2 + n_y*i1] = 1 + 0.5*c[i2 + n_y*i1];
			//__imag__ m[i2 + n_y*i1] = 0.0;
		}
	}

	gsl_rng_free(r);
	
	half_nx = (int) n_x/2;
	half_ny = (int) n_y/2;
	double n1real, n2real, creal;

	delta_kx = (2.0*M_PI)/(n_x*delta_x);
	delta_ky = (2.0*M_PI)/(n_y*delta_y);

	for(INDEX = 0; INDEX <=  10000; ++INDEX)
	{
		printf("%d\n", INDEX);
		for(i1=0; i1<n_x; ++i1)
		{
			for(i2=0; i2<n_y; ++i2)
			{
			  creal = (__real__ c[i2+n_y*i1]);
			  n1real = (__real__ n1[i2+n_y*i1]);
			  n2real = (__real__ n2[i2+n_y*i1]);

			  __real__ g[i2+n_y*i1] = 2*A*creal*(1-creal)*(1-creal) - 2*A*creal*creal*(1 - creal) + 2*creal*(0.25 + (-(n1real*n1real)/2 + (n1real*n1real*n1real*n1real)/2 - (n2real*n2real)/2 + (n2real*n2real*n2real*n2real)/2) + e*n1real*n1real*n2real*n2real);
			  __imag__ g[i2+n_y*i1] = 0.0;

			  __real__ h1[i2+n_y*i1] = (2 + creal*creal)*(-n1real + 2*n1real*n1real*n1real + 2*e*n1real*n2real);
			  __imag__ h1[i2+n_y*i1] = 0.0;

			  __real__ h2[i2+n_y*i1] = (2 + creal*creal)*(-n2real + 2*n2real*n2real*n2real + 2*e*n2real*n1real);
			  __imag__ h2[i2+n_y*i1] = 0.0;
			}
		}
		
		/** Let us take comp to the Fourier space **/

		fftw_execute_dft(planFc,c,c);
		fftw_execute_dft(planFg,g,g);
		fftw_execute_dft(planFn1,n1,n1);
		fftw_execute_dft(planFn2,n2,n2);
		fftw_execute_dft(planFh1,h1,h1);
		fftw_execute_dft(planFh2,h2,h2);
		fftw_execute_dft(planFctemp,ctemp,ctemp);
		fftw_execute_dft(planFn1temp,n1temp,n1temp);
		fftw_execute_dft(planFn2temp,n2temp,n2temp);

		/** Evolve composition **/

		for(i1=0; i1 < n_x; ++i1)
		{
			if(i1 < half_nx) 
				kx = i1*delta_kx;
			else 
				kx = (i1-n_x)*delta_kx;
			kx2 = kx*kx;

			for(i2=0; i2 < n_y; ++i2)
			{
				if(i2 < half_ny) 
					ky = i2*delta_ky;
				else 
					ky = (i2-n_y)*delta_ky;
				ky2 = ky*ky;
	
				k2 = kx2 + ky2;
				k4 = k2*k2;
				
				//coe1 = 1 + 2*k4*delta_t*(mAA*kA);
				//coe2 = 2*k4*delta_t*(mAB*kB);
				//coe3 = 2*k4*delta_t*(mAB*kA);
				//coe4 = 1 + 2*k4*delta_t*(mBB*kB);
				//det = coe1*coe4 - coe2*coe3;
				 
				ctemp[i2 + n_y*i1] = (c[i2+n_y*i1] - delta_t*M*k2*g[i2+n_y*i1])/(1+2*delta_t*M*kc*k4);
				
				n1temp[i2 + n_y*i1] = (n1[i2+n_y*i1] - delta_t*L*h1[i2+n_y*i1])/(1+2*delta_t*L*kn*k2);
				n2temp[i2 + n_y*i1] = (n2[i2+n_y*i1] - delta_t*L*h2[i2+n_y*i1])/(1+2*delta_t*L*kn*k2);
				
				c[i2 + n_y*i1] = ctemp[i2 + n_y*i1];
				n1[i2 + n_y*i1] = n1temp[i2 + n_y*i1];
				n2[i2 + n_y*i1] = n2temp[i2 + n_y*i1];
			}
		}

		/** Take composition back to real space **/

		fftw_execute_dft(planBc,c,c);
		fftw_execute_dft(planBg,g,g);
		fftw_execute_dft(planBn1,n1,n1);
		fftw_execute_dft(planBn2,n2,n2);
		fftw_execute_dft(planBh1,h1,h1);
		fftw_execute_dft(planBh2,h2,h2);
		fftw_execute_dft(planBctemp,ctemp,ctemp);
		fftw_execute_dft(planBn1temp,n1temp,n1temp);
		fftw_execute_dft(planBn2temp,n2temp,n2temp);
		
		for(i1 = 0; i1 < n_x; ++i1)
		{
			for(i2 = 0; i2 < n_y; ++i2)
			{
				c[i2+n_y*i1] = c[i2+n_y*i1]/(n_x*n_y); 
				n1[i2+n_y*i1] = n1[i2+n_y*i1]/(n_x*n_y); 
				n2[i2+n_y*i1] = n2[i2+n_y*i1]/(n_x*n_y); 
				__imag__ c[i2+n_y*i1] = 0.0;
				__imag__ n1[i2+n_y*i1] = 0.0;
				__imag__ n2[i2+n_y*i1] = 0.0;
			}
			
		}
		
		if(INDEX%TPRINT == 0)
		{
			sprintf(name,"./output/timeStep%d_c.dat", INDEX);
			fpo = fopen(name, "w");
			sprintf(name,"./output/timeStep%d_n1.dat", INDEX);
			fpon1 = fopen(name, "w");
			sprintf(name,"./output/timeStep%d_n2.dat", INDEX);
			fpon2 = fopen(name, "w");
						
			for(i1 = 0; i1 < n_x; ++i1)
			{
				for(i2 = 0; i2 < n_y; ++i2)
				{
				  fprintf(fpo, "%le ", __real__ c[i2+n_y*i1]);
				  fprintf(fpon1, "%le ", __real__ n1[i2+n_y*i1]);
				  fprintf(fpon2, "%le ", __real__ n2[i2+n_y*i1]);
				}
			        fprintf(fpo, "\n");
			        fprintf(fpon1, "\n");
			        fprintf(fpon2, "\n");
			}
			
			fprintf(gp, "plot \"%s\" matrix with image\n", name);
			fprintf(gp, "pause 0.05\n");
			fclose(fpo);
			fclose(fpon1);
			fclose(fpon2);
		}	
	}
	
	fftw_free(c);
	fftw_free(n1);
        fftw_free(n2);
	fftw_free(g);
	fftw_free(h1);
        fftw_free(h2);
	fftw_free(ctemp);
	fftw_free(n1temp);
        fftw_free(n2temp);
	fclose(gp);
	fftw_destroy_plan(planFc);
	fftw_destroy_plan(planFn1);
	fftw_destroy_plan(planFn2);
	fftw_destroy_plan(planFg);
	fftw_destroy_plan(planFh1);
	fftw_destroy_plan(planFh2);
	fftw_destroy_plan(planBc);
	fftw_destroy_plan(planBn1);
	fftw_destroy_plan(planBn2);
	fftw_destroy_plan(planBg);
	fftw_destroy_plan(planBh1);
	fftw_destroy_plan(planBh2);
}

