/*! \file steinhardt.h                                                                  
 *                                                                               
 * \brief Computes the Ql and Wl orientational order parameters
 * according to the definitions given in                        
 *      Steinhardt P. J., Nelson D.R. and Rochetti M (1983).    
 *      Phys. Rev. B, 28, 784                                   
 *      Bond-orientational order in liquids and glasses   
 * 		\url http://dx.doi.org/10.1103/PhysRevB.28.784         
 * This code takes advantage of the implementations of the  Legendre  
 * Polynomials and Wigner 3j Symbols provided in the GNU Scientific   
 * Library (GSL):                                                     
 * gsl_sf_legendre gsl_sf_coupling 
 * 
 * 
 * steindhardt.c is a small C library to calculate the Q_l and W_l 
 * orientational order parameters introduced by by Steindhardt and 
 * coworkers in Phys. Rev. B 28, 784 (1983) 
 * \url http://dx.doi.org/10.1103/PhysRevB.28.784
 * . It uses the GNU
 * Scientific Library to calculate the Legendre Polynomials 
 * (\url http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html) 
 * and Wigner 3j symbols 
 * (\url http://mathworld.wolfram.com/Wigner3j-Symbol.html) 
 * necessary to obtain Spherical Harmonics 
 * (\url http://mathworld.wolfram.com/SphericalHarmonic.html) and 
 * contract them in the W_l function. The code takes explicit advantage 
 * of the symmetry properties of the Spherical Harmonics to calculate 
 * only the ones with m>=0. The library provides 4 functions whose usage is 
 * explained below.
 */ 

/*! \fn int qlm (int l, int n, double rc, double *x, double *y, 
 * double *z, double *qlmRe, double *qlmIm):
 * This function calculates the quantities q_{lm} (0<=m<=l):
 * \f$ q_lm=\sum_{|r_{ij}|<rc} \mathcal{Y}_{lm}(\theta(r_{ij}),\phi(r_{ij})) \f$
 * The sum goes over all the possible bonds between the n particles 
 * whose coordinates are in the arrays x,y,z that are in magnitude 
 * smaller than rc. It will return the real and imaginary part of this 
 * quantities in the arrays qlmIm and qlmRe that have l+1 components. 
 * The function returns the 
 * value N of the number of bonds that fullfiled the condition 
 * \f$ |r_{ij}|<rc \f$.
 * Notice that in the literature the (\url http://dx.doi.org/10.1080/00268970210164428
 *  "C. Chakravarty Molecular Physics 2002, Vol 100, 
 * No 23, 37773780") the normalized quantities \f$ Q_{lm}=q_{lm}/N \f$ are the one 
 * that are actually calculated. The nomalization factor N is always
 * necessary for the calculation of the first order invariants
 * Q_l (which are calculated by the function Ql) but as long as N > 0
 * they are irrelevant for for the calculation of the second order
 * invariants W_l (which are calculated by the function Wl).
*/

/*! \fn double Ql (int l, int count, double *qlmRe, double *qlmIm)
 * This function calculates the *normalized* first order invariants:
 * \f$Q_l=\left(\frac{4\pi}{2l+1}\frac{1}{N} \sum_{m=-l}^{l}|q_{lm}| \right)^{1/2}
 * =\left(\frac{4\pi}{2l+1} \sum_{m=-l}^{l}|Q_{lm}| \right)^{1/2}
 * \f$ where the normalizing factor N, is given in the value of the 
 * variable count that is the integer value return by the function Qlm.
 * Also notice that because of the properties of the sphrerical 
 * harmonics only the Q_lm with positive or zero m are required to 
 * calculate the invariants.
 */ 

/*! \fn double wl (int l, double *qlmRe, double *qlmIm)
 * This function calculates the unnormalized second order invarinants
 * \f$w_l=\sum_{m_1+m_2+m_3} 
 * \left(\begin{array}{ccc} l & l & l \\ m_1 & m_2 & m_3 \end{array}\right)
 * q_{lm_1} q_{lm_2} q_{lm_3}\f$
 * where the the real and imaginary parts of the q_lm are provided
 * in the arrays qlmRe qlmIm that were previously calculated by the 
 * function qlm. The integer l indicate the value of the order
 * of the invariant and the dimensions of the arrays qlmRe and qlmIm
 * are l+1. The proper normalization of the sphrerical
 * harmonics is given \f$ W_l=w_l/q_{sum}^3 \f$ where \f$ q_{sum} \f$
 * is provided by the function qsum.
 */ 

/*! \fn double qsum (int l, double *qlmRe, double *qlmIm)
 * This function calculates \f$q_{sum}=\sqrt{\sum_{m=-l}^l |q_{lm}|^2}=
 * \frac{1}{N}\sqrt{\frac{2l+1}{4 \pi}} Q_l \f$ where again
 * N is the number of bonds which are returned by the function
 * qlm and Q_l is a first order invariant.
 * This quantity is used to normalize the second order invariants Wl
 */

int qlm (int l, int n, double rc, double *x, double *y, double *z,
	 double *qlmRe, double *qlmIm);
double qsum (int l, double *qlmRe, double *qlmIm);
double Ql (int l, int count, double *qlmRe, double *qlmIm);
double wl (int l, double *qlmRe, double *qlmIm);
