# Elliptic integrals, Jacobi's elliptic functions and Theta function.

The [Matlab](http://www.mathworks.com/)/[Octave](https://octave.org/) implementation of [Elliptic integrals of three types](https://en.wikipedia.org/wiki/Elliptic_integral), [Jacobi's elliptic functions](https://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi theta functions](https://en.wikipedia.org/wiki/Theta_function) of four types with their derivatives.

The main GOAL of the project is to provide the natural Matlab/Octave scripts WITHOUT external library calls like Maple and others. All scripts are developed to accept tensors as arguments and almost all of them have their complex versions. Performance and complete control on the execution are the main features.

Source code: [https://github.com/moiseevigor/elliptic](https://github.com/moiseevigor/elliptic)

## Installation

```matlab
git clone https://github.com/moiseevigor/elliptic.git
cd elliptic
setup    % adds src/ to the MATLAB/Octave path
```

# Contents of the package

  - [Elliptic Functions](elliptic#elliptic-functions)
    - [ELLIPJ: Jacobi's elliptic functions](elliptic#ellipj-jacobis-elliptic-functions)
    - [ELLIPJI: Jacobi's elliptic functions of the complex argument](elliptic#ellipji-jacobis-elliptic-functions-of-the-complex-argument)
    - [JACOBITHETAETA: Jacobi's Theta and Eta Functions](elliptic#jacobithetaeta-jacobis-theta-and-eta-functions)
    - [THETA: Theta Functions of Four Types](elliptic#theta-theta-functions-of-four-types)
    - [THETA_PRIME: Theta Functions and Their Derivatives](elliptic#theta_prime-theta-functions-and-their-derivatives)
  - [Elliptic Integrals](Elliptic-Integrals)
    - [ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function](elliptic#elliptic12)
    - [ELLIPTIC12I: Incomplete Elliptic Integrals of the complex argument](elliptic#elliptic12i)
    - [ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind](elliptic#elliptic3)
    - [ELLIPTIC123: Complete and Incomplete Elliptic Integrals](elliptic#elliptic123)
    - [INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind](elliptic#inverselliptic2)
  - [Weierstrass's elliptic functions (in development)](Weierstrass-Elliptic-Functions)
  - [Elliptic Related Functions](elliptic#elliptic-related-functions)
    - [AGM: Arithmetic Geometric Mean](elliptic#agm-arithmetic-geometric-mean)
    - [NOMEQ: The Value of Nome q = q(m)](elliptic#nomeq)
    - [INVERSENOMEQ: The Value of Nome m = m(q)](elliptic#inversenomeq)
  - [Jacobi Elliptic Functions (theory)](Jacobi-Elliptic-Functions)
  - [Theta Functions (theory)](Theta-Functions)
  - [Elliptic Curve](Elliptic-Curve)

# References

  1. [NIST Digital Library of Mathematical Functions](http://dlmf.nist.gov/)
  1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](http://www.nrbook.com/abramowitz_and_stegun/)" Dover Publications", 1965, Ch. 17.1 - 17.6.
  1. D. F. Lawden, "[Elliptic Functions and Applications](http://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989
  1. S. Zhang, J. Jin, "[Computation of Special Functions](http://jin.ece.uiuc.edu/specfunc.html)" (Wiley, 1996).
  1. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](http://nvl.nist.gov/pub/nistpubs/jres/107/5/j75car.pdf)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.
  1. N. H. Abel, "[Studies on Elliptic Functions](http://mathdl.maa.org/mathDL/46/?pa=content&sa=viewDocument&nodeId=1557)", english translation from french by Marcus Emmanuel Barnes. Original "Recherches sur les fonctions elliptiques", Journal fr die reine und angewandte Mathematik, Vol. 2, 1827. pp. 101-181.
  1. B. C. Berndt, H. H. Chan, S.-S. Huang, "[Incomplete Elliptic Integrals in Ramanujan's Lost Notebook](http://www.math.uiuc.edu/~berndt/articles/ellipticintegral.pdf)" in q-series from a Contemporary Perspective, M. E. H. Ismail and D. Stanton, eds., Amer. Math. Soc., 2000, pp. 79-126.
  1. J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)
