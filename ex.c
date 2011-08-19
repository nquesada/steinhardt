#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include"steinhardt.h"

int
main (int argv, char *argc[])
{
  int i, j;
  int tmp;
  int count;
  int n = atoi (argc[2]);
  int l = atoi (argc[3]);

  double rc = atof (argc[3]);
  double qss = 0;
  double *x, *y, *z;
  double *qlRe, *qlIm;
  double *qlval, *wlval;

  FILE *pf;

  x = (double *) malloc (n * sizeof (double));
  y = (double *) malloc (n * sizeof (double));
  z = (double *) malloc (n * sizeof (double));


  pf = fopen (argc[1], "r");
  if (pf != NULL && argv > 4)
    {

      for (i = 0; i < n; i++)
	{
	  tmp = fscanf (pf, "%lf %lf %lf", &x[i], &y[i], &z[i]);
	}
      fprintf (stdout, "\n%d\t", n);
      qlval = (double *) malloc ((argv - 4) * sizeof (double));
      wlval = (double *) malloc ((argv - 4) * sizeof (double));
      for (j = 4; j < argv; j++)
	{

	  l = atoi (argc[j]);
	  qlRe = (double *) malloc ((l + 1) * sizeof (double));
	  qlIm = (double *) malloc ((l + 1) * sizeof (double));

	  count = qlm (l, n, rc, x, y, z, qlRe, qlIm);
	  fprintf(stdout,"%d ",count);
	  qss = qsum (l, qlRe, qlIm);
	  qlval[j - 4] = Ql (l, count, qlRe, qlIm);
	  wlval[j - 4] = wl (l, qlRe, qlIm) / (qss * qss * qss);

	  free (qlRe);
	  free (qlIm);
	}
      for (j = 4; j < argv; j++)
	{
	  fprintf (stdout, "%.4lf\t", qlval[j - 4]);
	}
      for (j = 4; j < argv; j++)
	{
	  fprintf (stdout, "%.4lf\t", wlval[j - 4]);
	}
    }
  else
    {
      fprintf (stdout, "Could't find the file %s \n", argc[1]);
      exit (0);
    }
  free (x);
  free (y);
  free (z);
  fclose (pf);
  return 0;
}
