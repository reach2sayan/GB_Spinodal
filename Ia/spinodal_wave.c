#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

int main(void)
{

  FILE *fp10, *fp20, *fp50, *fpo;

  fp10 = fopen("./output/timeStep100_c.dat","r");
  fp20 = fopen("./output/timeStep200_c.dat","r");
  fp50 = fopen("./output/timeStep500_c.dat","r");
  fpo = fopen("spinodal_wave.dat", "w");

  int i1, i2;
  int n_x = 512, n_y = 256;
  double avg10, avg20, avg50;

  double *c10, *c20, *c50;
  
  c10 = malloc(n_x*n_y*sizeof(double));
  c20 = malloc(n_x*n_y*sizeof(double));
  c50 = malloc(n_x*n_y*sizeof(double));

  for(i1 = 0; i1 < n_x; i1++)
    {
      for(i2 = 0; i2 < n_y; i2++)
	{
	  fscanf(fp10, "%le ", &c10[i2+n_y*i1]);
	  fscanf(fp20, "%le ", &c20[i2+n_y*i1]);
	  fscanf(fp50, "%le ", &c50[i2+n_y*i1]);
	}
      fscanf(fp10, "\n");
      fscanf(fp20, "\n");
      fscanf(fp50, "\n");
    }
  fclose(fp10);
  fclose(fp20);
  fclose(fp50);

  avg10 = avg20 = avg50 = 0.0;

  for(i1 = (n_x/2)-50; i1 < (n_x/2)+50; i1++)
    {
      avg10 = avg20 = avg50 = 0.0;
      for(i2 = 0; i2 < n_y; i2++)
	{
	  avg10 = avg10 + c10[i2+n_y*i1];
	  avg20 = avg20 + c20[i2+n_y*i1];
	  avg50 = avg50 + c50[i2+n_y*i1];
	}
      fprintf(fpo,"%d %le %le %le\n", i1-(n_x)/2, (avg10)/n_y-0.5, (avg20)/n_y-0.5, (avg50)/n_y-0.5);
    }
}
