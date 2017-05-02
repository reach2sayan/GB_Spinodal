#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

int main(void)
{

  FILE *fp10, *fp20, *fp30, *fp40, *fp50, *fp60, *fp70, *fp80, *fp90, *fp100,*fpo;

  fp10 = fopen("./output/timeStep100_c.dat","r");
  fp30 = fopen("./output/timeStep300_c.dat","r");
  fp20 = fopen("./output/timeStep200_c.dat","r");
  fp40 = fopen("./output/timeStep400_c.dat","r");
  fp50 = fopen("./output/timeStep500_c.dat","r");
  fp60 = fopen("./output/timeStep600_c.dat","r");
  fp70 = fopen("./output/timeStep700_c.dat","r");
  fp80 = fopen("./output/timeStep800_c.dat","r");
  fp90 = fopen("./output/timeStep900_c.dat","r");
  fp100 = fopen("./output/timeStep1000_c.dat","r");
  fpo = fopen("spinodal_wave.dat", "w");

  int i1, i2;
  int n_x = 512, n_y = 256;
  double avg10, avg20, avg30, avg40, avg50, avg60, avg70, avg80, avg90, avg100;

  double *c10, *c20, *c30, *c40, *c50, *c60, *c70, *c80, *c90, *c100;
  
  c10 = malloc(n_x*n_y*sizeof(double));
  c20 = malloc(n_x*n_y*sizeof(double));
  c30 = malloc(n_x*n_y*sizeof(double));
  c40 = malloc(n_x*n_y*sizeof(double));
  c50 = malloc(n_x*n_y*sizeof(double));
  c60 = malloc(n_x*n_y*sizeof(double));
  c70 = malloc(n_x*n_y*sizeof(double));
  c80 = malloc(n_x*n_y*sizeof(double));
  c90 = malloc(n_x*n_y*sizeof(double));
  c100 = malloc(n_x*n_y*sizeof(double));

  for(i1 = 0; i1 < n_x; i1++)
    {
      for(i2 = 0; i2 < n_y; i2++)
	{
	  fscanf(fp10, "%le ", &c10[i2+n_y*i1]);
	  fscanf(fp20, "%le ", &c20[i2+n_y*i1]);
	  fscanf(fp30, "%le ", &c30[i2+n_y*i1]);
	  fscanf(fp40, "%le ", &c40[i2+n_y*i1]);
	  fscanf(fp50, "%le ", &c50[i2+n_y*i1]);
	  fscanf(fp60, "%le ", &c60[i2+n_y*i1]);
	  fscanf(fp70, "%le ", &c70[i2+n_y*i1]);
	  fscanf(fp80, "%le ", &c80[i2+n_y*i1]);
	  fscanf(fp90, "%le ", &c90[i2+n_y*i1]);
	  fscanf(fp100, "%le ", &c100[i2+n_y*i1]);
	}
      fscanf(fp10, "\n");
      fscanf(fp20, "\n");
      fscanf(fp30, "\n");
      fscanf(fp40, "\n");
      fscanf(fp50, "\n");
      fscanf(fp60, "\n");
      fscanf(fp70, "\n");
      fscanf(fp80, "\n");
      fscanf(fp90, "\n");
      fscanf(fp100, "\n");
    }
  fclose(fp10);
  fclose(fp20);
  fclose(fp30);
  fclose(fp40);
  fclose(fp50);
  fclose(fp60);
  fclose(fp70);
  fclose(fp80);
  fclose(fp90);
  fclose(fp100);

  avg10 = avg20 = avg30 = avg40 = avg50 = avg60 = avg70 = avg80 = avg90 = avg100 =0.0;

  for(i1 = (n_x/2)-50; i1 < (n_x/2)+50; i1++)
    {
      avg10 = avg20 = avg30 = avg40 = avg50 = avg60 = avg70 = avg80 = avg90 = avg100 =0.0;
      for(i2 = 0; i2 < n_y; i2++)
	{
	  avg10 = avg10 + c10[i2+n_y*i1];
	  avg20 = avg20 + c20[i2+n_y*i1];
	  avg30 = avg30 + c30[i2+n_y*i1];
	  avg40 = avg40 + c40[i2+n_y*i1];
	  avg50 = avg50 + c50[i2+n_y*i1];
	  avg60 = avg60 + c60[i2+n_y*i1];
	  avg70 = avg70 + c70[i2+n_y*i1];
	  avg80 = avg80 + c80[i2+n_y*i1];
	  avg90 = avg90 + c90[i2+n_y*i1];
	  avg100 = avg100 + c100[i2+n_y*i1];
	}
      fprintf(fpo,"%d %le %le %le %le %le %le %le %le %le %le\n", i1-(n_x)/2, (avg10)/n_y-0.5, (avg20)/n_y-0.5, (avg30)/n_y-0.5, (avg40)/n_y-0.5, (avg50)/n_y-0.5, (avg60)/n_y-0.5, (avg70)/n_y-0.5, (avg80)/n_y-0.5, (avg90)/n_y-0.5, (avg100)/n_y-0.5);
    }
}
