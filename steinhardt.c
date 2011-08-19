/*! \file steinhardt.c 
 * \brief Contains the source code. See documentation in steinhardt.h
 *
 * Copyright (C) 2010 Nicolas Quesada
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License ( GPL )
 * as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,
 *  USA.
 * http://www.gnu.org/licenses/gpl-2.0.html
 */

#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_sf_coupling.h>
#include<gsl/gsl_math.h>
#include<math.h>
#include"steinhardt.h"


int
qlm (int l, int n, double rc, double *x, double *y, double *z, double *qlmRe,
     double *qlmIm)
{
  //Variable Declaration
  int i, j, m;
  int count = 0;
  int lp1 = l + 1;

  double rc2 = rc * rc;
  double r2;
  double p;
  double xm, ym, zm, rm;
  double cost;
  double plm;
  for (i = 0; i < lp1; i++)
    {
      qlmRe[i] = 0.0;
      qlmIm[i] = 0.0;
    }

  for (i = 1; i < n; i++)
    {
      for (j = 0; j < i; j++)
	{
	  xm = x[i] - x[j];
	  ym = y[i] - y[j];
	  zm = z[i] - z[j];
	  r2 = xm * xm + ym * ym + zm * zm;
	  if (r2 < rc2)
	    {
	      count++;
	      rm = sqrt (r2);
	      cost = zm / rm;
	      p = atan2 (ym, xm);
	      plm = gsl_sf_legendre_sphPlm (l, 0, cost);
	      qlmRe[0] += plm;
	      for (m = 1; m <= l; m++)
		{
		  plm = gsl_sf_legendre_sphPlm (l, m, cost);
		  qlmRe[m] += plm * cos (m * p);
		  qlmIm[m] += plm * sin (m * p);
		}
	    }
	}
    }
  return count;
}



double
Ql (int l, int count, double *qlmRe, double *qlmIm)
{
  int m;
  double ql = 0;

  for (m = 1; m <= l; m++)
    {
      ql += qlmRe[m] * qlmRe[m] + qlmIm[m] * qlmIm[m];
    }
  ql = sqrt ((2 * ql + qlmRe[0] * qlmRe[0]) * 4 * M_PI / ((2.0 * l + 1)));

  ql = ql / count;

  return ql;

}


double
qsum (int l, double *qlmRe, double *qlmIm)
{
  int m;
  double ql = 0;

  for (m = 1; m <= l; m++)
    {
      ql += qlmRe[m] * qlmRe[m] + qlmIm[m] * qlmIm[m];
    }
  ql = sqrt (2 * ql + qlmRe[0] * qlmRe[0]);

  return ql;

}




double
wl (int l, double *qlmRe, double *qlmIm)
{
  int parity;
  int m1, m2, m3;
  int l2 = 2 * l;
  double x1, x2, x3, y1, y2, y3;
  double wl = 0;


  for (m1 = -l; m1 <= 0; m1++)
    {
      for (m2 = -l - m1; m2 <= l; m2++)
	{
	  m3 = -m1 - m2;

	  parity = 2 * GSL_IS_EVEN (m1) - 1;
	  x1 = parity * qlmRe[-m1];
	  y1 = -parity * qlmIm[-m1];


	  if (m2 < 0)
	    {
	      parity = 2 * GSL_IS_EVEN (m2) - 1;
	      x2 = parity * qlmRe[-m2];
	      y2 = -parity * qlmIm[-m2];
	    }
	  else
	    {
	      x2 = qlmRe[m2];
	      y2 = qlmIm[m2];
	    }
	  if (m3 < 0)
	    {
	      parity = 2 * GSL_IS_EVEN (m3) - 1;
	      x3 = parity * qlmRe[-m3];
	      y3 = -parity * qlmIm[-m3];
	    }
	  else
	    {
	      x3 = qlmRe[m3];
	      y3 = qlmIm[m3];
	    }
	  wl +=
	    gsl_sf_coupling_3j (l2, l2, l2, 2 * m1, 2 * m2,
				2 * m3) * (x1 * x2 * x3 - x1 * y2 * y3 -
					   x2 * y1 * y3 - x3 * y2 * y1);
	}
    }


  for (m1 = 1; m1 <= l; m1++)
    {
      for (m2 = -l; m2 <= l - m1; m2++)
	{
	  m3 = -m1 - m2;

	  x1 = qlmRe[m1];
	  y1 = qlmIm[m1];

	  if (m2 < 0)
	    {
	      parity = 2 * GSL_IS_EVEN (m2) - 1;
	      x2 = parity * qlmRe[-m2];
	      y2 = -parity * qlmIm[-m2];
	    }
	  else
	    {
	      x2 = qlmRe[m2];
	      y2 = qlmIm[m2];
	    }
	  if (m3 < 0)
	    {
	      parity = 2 * GSL_IS_EVEN (m3) - 1;
	      x3 = parity * qlmRe[-m3];
	      y3 = -parity * qlmIm[-m3];
	    }
	  else
	    {
	      x3 = qlmRe[m3];
	      y3 = qlmIm[m3];
	    }
	  wl +=
	    gsl_sf_coupling_3j (l2, l2, l2, 2 * m1, 2 * m2,
				2 * m3) * (x1 * x2 * x3 - x1 * y2 * y3 -
					   x2 * y1 * y3 - x3 * y2 * y1);
	}
    }
  return wl;
}
