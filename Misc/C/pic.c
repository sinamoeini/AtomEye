#include <stdio.h>
#include <math.h>
#define XMESH 500
#define YMESH 500

int main(int argc, char *argv[])
{
  FILE * IN;
  FILE * OUT;
  float  real[XMESH+1][YMESH+1];
  float  img [XMESH+1][YMESH+1];
  int i,j,k;

  if (argc!=2)
    {
      printf ("Input file name.\n");
      exit(0);
    }
  
  IN = fopen(argv[1],"r");
  OUT = fopen("a.out","w");

  for (i=0;i<=500;i+=10)
    {
      for (j=0;j<=500;j+=10)
	{
	  fscanf (IN,"%f %f\n",&real[i][j],&img[i][j]);
	  fprintf(OUT,"%f ",sqrt(real[i][j]*real[i][j]+img[i][j]
	     		   *img[i][j]));
	}
      fprintf(OUT,"\n");
    }
  
  fclose(IN);
  fclose(OUT);
  return(1);
}
    



