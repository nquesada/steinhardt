/*! \file \header steinhardt.h                                                                  
 *                                                                               
 * \brief Computes the Ql and Wl orientational order parameters.
 * 
 * \mainpage Orientational Order Parameters 
 * 
 * Computes the Ql and Wl orientational order parameters defined in:
 *      Steinhardt P. J., Nelson D.R. and Rochetti M (1983).    
 *      Phys. Rev. B, 28, 784                                   
 *      Bond-orientational order in liquids and glasses   
 * 		http://dx.doi.org/10.1103/PhysRevB.28.784         
 * This code takes advantage of the implementations of the  Legendre  
 * Polynomials and Wigner 3j Symbols provided in the GNU Scientific   
 * Library (GSL): gsl_sf_legendre and gsl_sf_coupling .
 * To see an example of the use of the library see the documentation of
 * example.c .
 * 
 * 
 * steindhardt is a small C library to calculate the \f$ Q_l \f$ and \f$
 * W_l \f$ orientational order parameters introduced by 
 * P. J. Steindhardt and coworkers in Phys. Rev. B 28, 784 (1983) 
 * http://dx.doi.org/10.1103/PhysRevB.28.784 . It uses the GNU
 * Scientific Library to calculate the Legendre Polynomials 
 * (http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html) 
 * and Wigner 3j symbols 
 * (http://mathworld.wolfram.com/Wigner3j-Symbol.html) 
 * necessary to obtain Spherical Harmonics 
 * (http://mathworld.wolfram.com/SphericalHarmonic.html) and 
 * contract them in rotationally invariant combinations. 
 * The code takes explicit advantage 
 * of the symmetry properties of the Spherical Harmonics to calculate 
 * only the ones with m>=0. The library provides 4 functions whose usage
 * is explained below.
 */ 

/*! \fn int qlm (int l, int n, double rc, double *x, double *y, 
 * double *z, double *qlmRe, double *qlmIm):
 * This function calculates the quantities q_lm (0<=m<=l):
 * \f$ q_{lm}=\sum_{|r_{ij}|<r_c} 
 * 		\mathcal{Y}_{lm}(\theta(r_{ij}),\phi(r_{ij})) \f$
 * The sum goes over all the possible bonds between the n particles 
 * whose coordinates are in the arrays x,y,z that are in magnitude 
 * smaller than r_c. \n
 * 
 * Input:	int l, integer of the order parameter. \n
 * 			int n, Number of particles of the cluster. \n
 * 			double rc, cutoff radius. Only bonds of length small that rc
 * 				will be considered. \n
 * 			double *x, *y, *z, the x, y and z coordinates of the atoms
 * 				of the cluster. Each of them is an array of dimension n.
 * 				\n
 * 			double *qlmRe, *qlIm these arrays will be used to return 
 * 				the values of the quantities qlm as defined in the above
 * 				equation. \n
 * 
 * Output:	The function returns the value of the number of bonds that
 * 			satisfied the cutoff condition |r_{ij}|<r_c .\n 
 * 			It will also return the real and imaginary part of this 
 * 			quantities (qlm) in the arrays qlmRe and qlIm that have 
 * 			l+1 components corresponding to the values of m from 0 to l. 
 *	 		Notice that in the literature C. Chakravarty Molecular 
 * 			Physics 2002, Vol 100, No 23, 37773780 
 * 			http://dx.doi.org/10.1080/00268970210164428 the normalized 
 * 			quantities \f$ Q_{lm}=q_{lm}/N \f$ are the ones 
 * 			that are actually calculated. The normalization factor N is 
 * 			always necessary for the calculation of the second order
 * 			invariants Q_l (which are calculated by the function Ql) 
 * 			but as long as N > 0 they are irrelevant for the calculation
 * 			of the third order invariants W_l 
 * 			(which are calculated by the function wl and qsum).
*/

/*! \fn double Ql (int l, int count, double *qlmRe, double *qlmIm)
 * This function calculates the *normalized* first order invariants:
 * \f$Q_l=\frac{1}{N}\left(\frac{4\pi}{2l+1} 
 * \sum_{m=-l}^{l}|q_{lm}|^2 \right)^{1/2}
 * =\left(\frac{4\pi}{2l+1} \sum_{m=-l}^{l}|Q_{lm}|^2 \right)^{1/2}
 * \f$ where the normalizing factor N, is given in the value of the 
 * variable count that is the integer value returned by the function 
 * qlm.
 * 
 * Input:	int l, integer of the order parameter. \n
 * 			int count, Number of bonds considered in the calculation of
 * 				the qlm. It is the value returned by the function qlm.\n
 * 			double *qlmRe, *qlmIm Real and imaginary parts of the l+1
 * 				numbers qlm that were calculated by the function qlm.\n
 * 
 * Output: 	Normalized value of the second order parameter (also 
 * 			called second order invariant) Ql.
 * 
 * Notice that because of the properties of the spherical 
 * harmonics only the qlm with positive or zero m are required to 
 * calculate the invariants.
 */ 

/*! \fn double wl (int l, double *qlmRe, double *qlmIm)
 * This function calculates the unnormalized second order invarinants
 * \f$w_l=\sum_{m_1+m_2+m_3=0} 
 * \left(\begin{array}{ccc} l & l & l \\ m_1 & m_2 & m_3 
 * \end{array}\right) q_{lm_1} q_{lm_2} q_{lm_3}\f$. \n
 * 
 * Input:	int l, integer of the order parameter..
 * 			double *qlmRe, *qlmIm Real and imaginary parts of the l+1
 * 			numbers qlm that were calculated by the function qlm.\n
 * 
 * Output: 	Unnormalized value of the third order parameter (also 
 * 			called third order invariant) wl.
 *			The proper normalization of the Wl is given by 
 *			\f$ W_l=w_l/q_{sum}^3 \f$ where \f$ q_{sum} \f$ is provided 
 *			by the function qsum.
 */ 

/*! \fn double qsum (int l, double *qlmRe, double *qlmIm)
 * This function calculates \f$q_{sum}=\sqrt{\sum_{m=-l}^l |q_{lm}|^2}
 * \f$. This quantity is used to normalize the second order invariants 
 * Wl, \f$ W_l=w_l/q_{sum}^3 \f$ where \f$ q_{sum} \f$. \n
 * Input:	int l, integer of the order parameter..
 * 			double *qlmRe, *qlmIm Real and imaginary parts of the l+1
 * 			numbers qlm that were calculated by the function qlm.\n
 * Output:	qsum as given in the above equation.
 */

int qlm (int l, int n, double rc, double *x, double *y, double *z,
	 double *qlmRe, double *qlmIm);
double qsum (int l, double *qlmRe, double *qlmIm);
double Ql (int l, int count, double *qlmRe, double *qlmIm);
double wl (int l, double *qlmRe, double *qlmIm);
