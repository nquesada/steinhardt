/*! \file example.c
 * \brief Provides and example of the usage of the library steinhardt.
 * 
 * 
 * Provides and example of the usage of the library steinhardt.
 * It also can be used to reproduce the results presented in 
 * ( C. Chakravarty Molecular Physics 2002, Vol 100, No 23, 37773780 
 * 	http://dx.doi.org/10.1080/00268970210164428) with the help of the
 * scripts getdata.sh and runexample.sh .
 * 
 * \fn int main (int argv, char *argc[])
 * The function reads from the standard input a coordinates file and 
 * calculates the order parameters Q_l and W_l.
 * To compile the program type in a terminal:\n
 *  make example.out
 * 
 * To execute the program the following arguments should be passed at
 * execution in the command line. One after the other in the following
 * order:\n
 * 		N, the number of atoms.\n
 *		file.xyz, name of the file containing the geometry. The first N 
 * 			lines of the file will be read.\n
 * 		rc, cutoff radius only bonds of length less than rc will be
 * 			considered in the calculation of the order parameters.\n
 * 		l1, l2, ... The following arguments should be integers 
 * 		specifying the order of the order parameters that the user
 * 		wishes to calculate.\n
 * 
 * 		An example:\n
 * 		./example.out 38 38DD.xyz 1.25 2 4 6 8 10 
 * 		\n
 * 		This specifies a cluster of 38 atoms with geometry in the 
 * 		file 38DD.xyz . The cutoff distance will be 1.25 and the order
 * 		parameters that will be calculated are l=2,4,6,8,10.
 * 		The geometry file should contain in 3 columns the x, y and z
 * 		components of the geometry and NOTHING else. This is the same
 * 		format used in the  The Cambridge Cluster Database, 
 * 		http://www-wales.ch.cam.ac.uk/CCD.html for Lennard-Jones
 * 		clusters.\n
 * 
 * The file will return separated by tabs the following:\n
 * Number of atoms, the values of the Q_l and then the values W_l
 * in the same order of the l values given as inputs (in the example
 * this would be 2,4,6,8,10).\n
 * 
 * Finally, the scrips getdata.sh and runexample.sh can be used to
 * fetch datafiles from the CCD and then use the program to generate 
 * the tables presented in C. Chakravarty Molecular Physics 2002, 
 * Vol 100, No 23, 37773780 http://dx.doi.org/10.1080/00268970210164428:
 * \n
 * Compile the executable:\n
 * 		make example.out\n
 * Fetch the data file:\n
 * 		bash getdata.sh\n
 * Generate the table:\n
 * 		bash runexample.sh\n
 */ 

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
	  qss = qsum (l, qlRe, qlIm);
	  if (qss > 1e-3)
	    {
	      qlval[j - 4] = Ql (l, count, qlRe, qlIm);
	      wlval[j - 4] = wl (l, qlRe, qlIm) / (qss * qss * qss);
	    }
	  else
	    {
	      qlval[j - 4] = 0.0;
	      wlval[j - 4] = 0.0;
	    }

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
