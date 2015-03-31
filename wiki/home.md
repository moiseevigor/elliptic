# Elliptic integrals, Jacobi's elliptic functions and Theta function.

<wiki:comment>
<wiki:gadget url="http://elliptic.googlecode.com/svn/trunk/gadgets/elliptic-gadget-1.xml" border="0" width="728" height="15" />
<wiki:gadget url="http://elliptic.googlecode.com/svn/trunk/gadgets/elliptic-gadget-3.xml" border="0" width="728" height="90" />
</wiki:comment>

The [Matlab](http://www.mathworks.com/) script implementations of [Elliptic integrals of three types](http://code.google.com/p/elliptic/wiki/EllipticIntegrals), [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi theta functions](http://en.wikipedia.org/wiki/Theta_function) of four types.

The main GOAL of the project is to provide the natural Matlab scripts WITHOUT external library calls like Maple and others. All scripts are developed to accept tensors as arguments and almost all of them have their complex versions. Performance and complete control on the execution are the main features.

To download please go to [Routines list](https://elliptic.googlecode.com/svn/trunk/).

Or checkout via [Subversion](http://subversion.tigris.org/) from [Checkout page](https://code.google.com/p/elliptic/source/checkout).

<wiki:comment>
[http://elliptic.googlecode.com/svn/trunk/imgs/elliptic.png](http://code.google.com/p/elliptic/wiki/elliptic)

**[Elliptic integrals](http://code.google.com/p/elliptic/wiki/EllipticIntegrals)** originally arose in connection with the problem of giving the arc length of an ellipse. They were first studied by Giulio Fagnano and Leonhard Euler. Modern mathematics defines an elliptic integral as any function f which can be expressed in the form
```
f(x) = Integral(R(t,P(t), c, x)dt,
```
where `R` is a rational function of its two arguments, `P` is the square root of a polynomial of degree `3` or `4` with no repeated roots, and `c` is a constant.
In general, elliptic integrals cannot be expressed in terms of elementary functions. Exceptions to this general rule are when `P` has repeated roots, or when `R(x,y)` contains no odd powers of `y`. However, with the appropriate reduction formula, every elliptic integral can be brought into a form that involves integrals over rational functions and the three canonical forms (i.e. the elliptic integrals of the first, second and third kind).

**[The Jacobi's elliptic functions](http://code.google.com/p/elliptic/wiki/elliptic#Elliptic_Functions)** are a set of basic elliptic functions, and auxiliary theta functions, that have historical importance with also many features that show up important structure, and have direct relevance to some applications (e.g. the equation of the pendulum). They also have useful analogies to the functions of trigonometry, as indicated by the matching notation `SN` for `SIN`. They are not the simplest way to develop a general theory, as now seen: that can be said for the Weierstrass elliptic functions. They are not, however, outmoded. They were introduced by Carl Gustav Jakob Jacobi, around 1830.

**[Theta functions](http://code.google.com/p/elliptic/wiki/elliptic#THETA:_Theta_Functions_of_Four_Types)** are special functions of several complex variables. They are important in several areas, including the theories of abelian varieties and moduli spaces, and of quadratic forms. They have also been applied to soliton theory. When generalized to a Grassmann algebra, they also appear in quantum field theory, specifically string theory and D-branes.

**[Weierstrass's elliptic functions](http://code.google.com/p/elliptic/wiki/elliptic#Weierstrass's_elliptic_functions_(in_development))** are elliptic functions that take a particularly simple form (cf Jacobi's elliptic functions); they are named for Karl Weierstrass. This class of functions are also referred to as p-functions and generally written using the symbol ℘ (a stylised letter p called Weierstrass p).

The Weierstrass elliptic function can be defined in three closely related ways, each of which possesses certain advantages. One is as a function of a complex variable z and a lattice Λ in the complex plane. Another is in terms of z and two complex numbers ω1 and ω2 defining a pair of generators, or periods, for the lattice. The third is in terms z and of a modulus τ in the upper half-plane. This is related to the previous definition by τ = ω2 / ω1, which by the conventional choice on the pair of periods is in the upper half-plane. Using this approach, for fixed z the Weierstrass functions become modular functions of τ.

<table border="0">
<tr>
<td valign="top" align="left" width="80px">
<a rel="nofollow" href="https://www.paypal.com/cgi-bin/webscr?cmd=_xclick&business=moiseev.igor@gmail.com&item_name=elliptic&no_shipping=0&no_note=1&tax=0&currency_code=EUR&lc=IT&bn=PP-DonationsBF&charset=UTF-8">
<img width="62" height="31" src="https://www.paypal.com/en_US/i/btn/x-click-but04.gif" alt="PayPal — donate"/>
</a>
</td>
<td valign="top" >
<div>Any donation will be appreciated!</div>
Thank you.
</td>
</tr>
</table>

</wiki:comment>

# Contents of the package

  - [Elliptic Functions](http://code.google.com/p/elliptic/wiki/elliptic#Elliptic_Functions)
    - [ELLIPJ: Jacobi's elliptic functions](http://code.google.com/p/elliptic/wiki/elliptic#ELLIPJ:_Jacobi's_elliptic_functions)
    - [ELLIPJI: Jacobi's elliptic functions of the complex argument](http://code.google.com/p/elliptic/wiki/elliptic?#ELLIPJI:_Jacobi's_elliptic_functions_of_the_complex_argumen)
    - [JACOBITHETAETA: Jacobi's Theta and Eta Functions](http://code.google.com/p/elliptic/wiki/elliptic#JACOBITHETAETA:_Jacobi's_Theta_and_Eta_Functions)
    - [THETA: Theta Functions of Four Types](http://code.google.com/p/elliptic/wiki/elliptic#THETA:_Theta_Functions_of_Four_Types)
  - [Elliptic Integrals](http://code.google.com/p/elliptic/wiki/elliptic#Elliptic_Integrals)
    - [ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function](http://code.google.com/p/elliptic/wiki/elliptic#ELLIPTIC12:_Incomplete_Elliptic_Integrals_of_the_First,_Second_K)
    - [ELLIPTIC12I: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function of the complex argument](http://code.google.com/p/elliptic/wiki/elliptic#ELLIPTIC12I:_Incomplete_Elliptic_Integrals_of_the_First,_Second)
    - [ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind](http://code.google.com/p/elliptic/wiki/elliptic#ELLIPTIC3:_Incomplete_Elliptic_Integral_of_the_Third_Kind)
    - [ELLIPTIC123: Complete and Incomplete Elliptic Integrals of the First, Second, and Third Kind](http://code.google.com/p/elliptic/wiki/elliptic#ELLIPTIC123:_Complete_and_Incomplete_Elliptic_Integrals_of_the_F)
    - <font color="red"></font> [INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind](http://code.google.com/p/elliptic/wiki/elliptic#INVERSELLIPTIC2:_INVERSE_Incomplete_Elliptic_Integrals_of_the_Se)

  - [Weierstrass's elliptic functions (in development)](http://code.google.com/p/elliptic/wiki/elliptic#Weierstrass's_elliptic_functions_(in_development))
  - [Elliptic Related Functions](http://code.google.com/p/elliptic/wiki/elliptic#Elliptic_Related_Functions)
    - [AGM: Artihmetic Geometric Mean](http://code.google.com/p/elliptic/wiki/elliptic#AGM:_Artihmetic_Geometric_Mean)
    - [NOMEQ: The Value of Nome q = q(m)](http://code.google.com/p/elliptic/wiki/elliptic#NOMEQ:_The_Value_of_Nome_q_=_q(m))
    - [INVERSENOMEQ: The Value of Nome m = m(q)](http://code.google.com/p/elliptic/wiki/elliptic#INVERSENOMEQ:_The_Value_of_Nome_m_=_m(q))
  - [References](http://code.google.com/p/elliptic/wiki/elliptic#References)

<wiki:comment>
# Contents of package

<wiki:toc max_depth="3" />

## Elliptic Functions

### ELLIPJ: Jacobi's elliptic functions

[ELLIPJ](http://code.google.com/p/elliptic/source/browse/trunk/ellipj.m) evaluates the [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi's amplitude](http://mathworld.wolfram.com/JacobiAmplitude.html).

`[= ELLIPJ(U,M)` returns the values of the Jacobi elliptic functions `SN`, `CN`, `DN` and `AM` evaluated for corresponding elements of argument U and parameter M.  The arrays U and M must be of the same size (or either can be scalar).  As currently implemented, M is limited to `0 <= M <= 1`. 

*General definition:*
```
u = Integral(1/sqrt(1-m^2*sin(theta)^2), 0, phi);
Sn(u) = sin(phi);
Cn(u) = cos(phi);
Dn(u) = sqrt(1-m^2*sin(phi)^2);
```

_Depends on_  `AGM`, `ELLIPKE`.<br>
_Used by_ `THETA`.<br>
_See also_ `ELLIPKE`.


### ELLIPJI: Jacobi's elliptic functions of the complex argument

[http://code.google.com/p/elliptic/source/browse/trunk/ellipji.m ELLIPJI](Sn,Cn,Dn,Am]) evaluates the Jacobi elliptic functions of complex phase `U`.

`[= ELLIPJ(U,M)` returns the values of the Jacobi elliptic functions `SNI`, `CNI` and `DNI` evaluated for corresponding  elements of argument `U` and parameter `M`. The arrays `U` and `M` must  be of the same size (or either can be scalar).  As currently implemented, `M` is real and limited to `0 <= M <= 1`. 


*Example:
```
[phi1,phi2](Sni,Cni,Dni]) = meshgrid(-pi:3/20:pi, -pi:3/20:pi);
phi = phi1 + phi2*i;
[= ellipji(phi, 0.99);
```

_Depends on_ `AGM`, `ELLIPJ`, `ELLIPKE`<br>
_See also_ `ELLIPTIC12`, `ELLIPTIC12I`

### JACOBITHETAETA: Jacobi's Theta and Eta Functions

[http://code.google.com/p/elliptic/source/browse/trunk/jacobiThetaEta.m JACOBITHETAETA](Sni,Cni,Dni]) evaluates Jacobi's theta and eta functions.

`[H](Th,) = JACOBITHETAETA(U,M)` returns the values of the Jacobi's theta and eta elliptic functions `TH` and `H` evaluated for corresponding elements of argument `U` and parameter `M`.  The arrays `U` and `M` must be the same size (or either can be scalar).  As currently implemented, `M` is real and limited to `0 <= M <= 1`. 

*Example:
```
[= meshgrid(0:5:90, 0:2:90);                  
[Th, H](phi,alpha]) = jacobiThetaEta(pi/180*phi, sin(pi/180*alpha).^2);  
```

*Depends on* `AGM`, `ELLIPJ`, `ELLIPKE`<br>
*See also* `ELLIPTIC12`, `ELLIPTIC12I`, `THETA` 
 
### THETA: Theta Functions of Four Types

[THETA](http://code.google.com/p/elliptic/source/browse/trunk/theta.m) evaluates theta functions of four types.

`Th = THETA(TYPE,V,M)` returns values of theta functions
evaluated for corresponding values of argument `V` and parameter `M`. `TYPE` is a type of the theta function, there are four numbered types. The arrays `V` and `M` must be the same size (or either can be scalar). As currently implemented, `M` is limited to `0 <= M <= 1`. 

`Th = THETA(TYPE,V,M,TOL)` computes the theta and eta elliptic functions to the accuracy `TOL` instead of the default `TOL = EPS`.  

The parameter `M` is related to the nome `Q` as `Q = exp(-pi*K(1-M)/K(M))`. Some definitions of the Jacobi's elliptic functions use the modulus `k` instead of the parameter `m`.  They are related by `m = k^2`.

*Example:
```
[= meshgrid(0:5:90, 0:2:90);                  
Th1 = theta(1, pi/180*phi, sin(pi/180*alpha).^2);  
Th2 = theta(2, pi/180*phi, sin(pi/180*alpha).^2);  
Th3 = theta(3, pi/180*phi, sin(pi/180*alpha).^2);  
Th4 = theta(4, pi/180*phi, sin(pi/180*alpha).^2);  
```

_Depends on_ `AGM`, `ELLIPJ`, `ELLIPKE`, `JACOBITHETAETA`<br>
_See also_ `ELLIPTIC12`, `ELLIPTIC12I`

## Elliptic Intergrals

### ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function

[http://code.google.com/p/elliptic/source/browse/trunk/elliptic12.m ELLIPTIC12](phi,alpha]) evaluates the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function.

`[= ELLIPTIC12(U,M,TOL)` uses the method of the Arithmetic-Geometric Mean and Descending Landen Transformation described in [http://code.google.com/p/elliptic/#References 1](F,E,Z]) Ch. 17.6, to determine the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function (see [1](http://code.google.com/p/elliptic/#References), [2](http://code.google.com/p/elliptic/#References)).

**General definition:**
```
F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
```


Tables generating code (see [1](http://code.google.com/p/elliptic/#References), pp. 613-621):
```
[= meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees
[F,E,Z](phi,alpha]) = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals
```

*Depends on* `AGM`<br>
*See also* `ELLIPKE`, `ELLIPJ`, `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

### ELLIPTIC12I: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function of the complex argument

[ELLIPTIC12i](http://code.google.com/p/elliptic/source/browse/trunk/elliptic12i.m) evaluates the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function for the complex value of phase `U`. Parameter `M` must be in the range `0 <= M <= 1`. 

`[= ELLIPTIC12i(U,M,TOL)` where `U` is a complex phase in radians, `M` is the real parameter and `TOL` is the tolerance (optional). Default value for the tolerance is `eps = 2.220e-16`.

`ELLIPTIC12i` uses the function `ELLIPTIC12` to evaluate the values of corresponding integrals.

*Example:
```
[phi1,phi2](Fi,Ei,Zi]) = meshgrid(-2*pi:3/20:2*pi, -2*pi:3/20:2*pi);
phi = phi1 + phi2*i;
[= elliptic12i(phi, 0.5);
```

_Depends on_ `ELLIPTIC12`, `AGM`<br>
_See also_ `ELLIPKE`, `ELLIPJ`, `ELLIPTIC3`, `THETA`.
 
### ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind

[http://code.google.com/p/elliptic/source/browse/trunk/elliptic3.m ELLIPTIC3](Fi,Ei,Zi]) evaluates incomplete elliptic integral of the third kind `Pi = ELLIPTIC3(U,M,C)` where `U` is a phase in radians, `0 < M < 1` is the module and `0 < C < 1` is a parameter. 

`ELLIPTIC3` uses Gauss-Legendre 10 points quadrature template described in `[to determine the value of the Incomplete Elliptic Integral of the Third Kind (see `[1, 2](3]`)`).

*General definition:
```
Pi(u,m,c) = int(1/((1 - c*sin(t)^2)*sqrt(1 - m*sin(t)^2)), t=0..u)
```

Tables generating code ([pp. 625-626):
```
[phi,alpha,c](1],) = meshgrid(0:15:90, 0:15:90, 0:0.1:1);
Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
```

## Weierstrass's elliptic functions (in development)

<font color="red">!IN DEVELOPMENT, help needed!</font>

## Elliptic Related Functions

### AGM: Artihmetic Geometric Mean

[AGM](http://code.google.com/p/elliptic/source/browse/trunk/agm.m) calculates the [Artihmetic Geometric Mean](http://en.wikipedia.org/wiki/Arithmetic-geometric_mean) of `A` and `B` (see [1](http://code.google.com/p/elliptic/#References)). 

`[= AGM(A0,B0,C0,TOL)` carry out the process of the arithmetic geometric mean, starting with a given positive numbers triple `(A0, B0, C0)` and returns in 
`(A, B, C)` the generated sequence. `N` is a number of steps (returns in the value`uint32`).

The general scheme of the procedure:
```
A(i) = 1/2*( A(i-1)+B(i-1) );     A(0) = A0;
B(i) = sqrt( A(i-1)*B(i-1) );     B(0) = B0;
C(i) = 1/2*( A(i-1)+B(i-1) );     C(0) = C0;
```
Stop at the `N`-th step when `A(N) = B(N)`, i.e., when `C(N) = 0`. 

_Used by_  `ELLIPJ` and `ELLIPTIC12`.<br>_See also_ `ELLIPKE`, `ELLIPTIC3`, `THETA`.

### NOMEQ: The Value of Nome `q = q(m)`

[http://code.google.com/p/elliptic/source/browse/trunk/nomeq.m NOMEQ](A,B,C,N]) gives the value of Nome `q = q(m)`.

Nome `Q = nomeq(M,TOL)`, where `0<=M<=1` is the module and `TOL` is the tolerance (optional). Default value for the tolerance is `eps = 2.220e-16`.

*Used by*  `ELLIPJ`.<br>
*Depends on* `ELLIPKE`<br>
*See also* `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

### INVERSENOMEQ: The Value of Nome `m = m(q)`

[INVERSENOMEQ](http://code.google.com/p/elliptic/source/browse/trunk/inversenomeq.m) gives the value of Nome `m = m(q)`.

`M = inversenomeq(q)`, where `Q` is the Nome of q-series.

**WARNING**. The function `INVERSENOMEQ` does not return correct values of `M` for `Q > 0.6`, because of computer precision limitation. The function `NomeQ(m)` has an essential singularity at `M = 1`, so it cannot be inverted at this point and actually it is very hard to find and inverse in the neigborhood also.

More preciesly:
```
nomeq(1) = 1
nomeq(1-eps) = 0.77548641878026
```

*Example:
```
nomeq(inversenomeq([0.3 0.4 0.5 0.6 0.7 0.8](0.001)))
```

*Used by*  `ELLIPJ`.<br>
*Depends on* `ELLIPKE`<br>
*See also* `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

# References

  1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](http://www.math.ucla.edu/~cbm/aands/)" Dover Publications", 1965, Ch. 17.1 - 17.6.
  1. D. F. Lawden, "[Elliptic Functions and Applications](http://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989
  1. S. Zhang, J. Jin, "[Computation of Special Functions](http://jin.ece.uiuc.edu/specfunc.html)" (Wiley, 1996).
  1. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](http://nvl.nist.gov/pub/nistpubs/jres/107/5/j75car.pdf)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.
  1. N. H. Abel, "[Studies on Elliptic Functions](http://mathdl.maa.org/mathDL/46/?pa=content&sa=viewDocument&nodeId=1557)", english translation from french by Marcus Emmanuel Barnes. Original "Recherches sur les fonctions elliptiques", Journal fr die reine und angewandte Mathematik, Vol. 2, 1827. pp. 101-181.
  1. B. C. Berndt, H. H. Chan, S.-S. Huang, "[Incomplete Elliptic Integrals in Ramanujan's Lost Notebook](http://www.math.uiuc.edu/~berndt/articles/ellipticintegral.pdf)" in q-series from a Contemporary Perspective, M. E. H. Ismail and D. Stanton, eds., Amer. Math. Soc., 2000, pp. 79-126.
</wiki:comment>

<wiki:comment>

# External Links

  - [Elliptic integrals](http://en.wikipedia.org/wiki/Elliptic_integral) and [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi's_elliptic_functions) from Wikipedia.
  - [ALGLIB](http://www.alglib.net/) - cross-platform numerical analysis and data processing library.
  - [Graphs by Kojocho2](http://www.youtube.com/user/kojocho2), Graph of complex functions and complex analysis.

# Some interesting applications

### Fingerprints In The Sky Explained By 'Beautiful Mathematics'.  [ScienceDaily (Nov. 23, 2004)](http://www.sciencedaily.com/releases/2004/11/041123112414.htm)

  *This is beautiful mathematics in the sky. Using elliptic integrals, we’ve been able to replace pages and pages of formulae with one very simple solution that predicts the pattern extremely well.* 

### Mathematicians Solve 'Trillion Triangle' Problem.  [ScienceDaily (Sep. 22, 2009)](http://www.sciencedaily.com/releases/2009/09/090922095651.htm)

  *According to Brian Conrey, Director of the American Institute of Mathematics, "Old problems like this may seem obscure, but they generate a lot of interesting and useful research as people develop new ways to attack them."*

### Elliptic-Function Filter Properties. [Connexions](http://cnx.org/content/m16925/latest/)

  *...In order to develop analytical expressions for equal-ripple rational functions, an interesting class of transcendental functions, called the Jacobian elliptic functions, is outlined.*
</wiki:comment>

<wiki:comment>
<br>
<br>
<wiki:gadget url="http://elliptic.googlecode.com/svn/trunk/gadgets/elliptic-gadget-1.xml" border="0" width="728" height="15" />
</wiki:comment>