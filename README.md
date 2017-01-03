# Computes the Ql and Wl orientational order parameters.

Computes the Ql and Wl orientational order parameters defined in:


       Steinhardt P. J., Nelson D.R. and Rochetti M.    
       Phys. Rev. B, 28.2 (1983): 784                                   
       Bond-orientational order in liquids and glasses   
       http://dx.doi.org/10.1103/PhysRevB.28.784


This code takes advantage of the implementations of the  [Legendre Polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials) and [Wigner 3j Symbols](https://en.wikipedia.org/wiki/Legendre_polynomials) provided in the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/): [`gsl_sf_legendre`](https://www.gnu.org/software/gsl/manual/html_node/Legendre-Functions-and-Spherical-Harmonics.html#Legendre-Functions-and-Spherical-Harmonics) and [`gsl_sf_coupling`](https://www.gnu.org/software/gsl/manual/html_node/Coupling-Coefficients.html#Coupling-Coefficients) .

A good reference on orientational order parameters is:


       Chakravarty C.
       Mol. Phys. 100.23 (2002): 3777-3780
       Bond orientational order in atomic clusters
       http://dx.doi.org/10.1080/00268970210164428

To see the documentation of the code see the source code of each file.
The documentation can be automatically extracted using Doxygen.

To compile an example type:
```bash
make example.out
```
To fetch the geometry files from the Cambridge Cluster Database type:
```bash
make data
```
To use these files to generate the tables cited in the documentation type:
```bash
make run
```
