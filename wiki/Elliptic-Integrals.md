# Elliptic integrals, numerical methods, Matlab examples of usage and problem solving.

<wiki:comment>
<wiki:gadget url="http://elliptic.googlecode.com/svn/trunk/gadgets/elliptic-gadget-1.xml" border="0" width="728" height="15" />
</wiki:comment>

<wiki:toc max_depth="3" />

# Elliptic Integrals

If `R(x,y)` is a rational function of `x` and `y`, where `y`^2^ is equal to a cubic or quadratic polynomial in `x`, the integral

```
f(x) = Integral(R(x,y))dx
```

is called an **elliptic integral**.
Elliptic integrals originally arose in the connection with the problem of giving the arc length of an ellipse. They were first studied by Giulio Fagnano and Leonhard Euler.

In general, elliptic integrals cannot be expressed in terms of elementary functions. Exceptions to this general rule are when `P` has repeated roots, or when `R(x,y)` contains no odd powers of `y`. However, with the appropriate reduction formula, every elliptic integral can be brought into a form that involves integrals over rational functions and the three canonical forms (i.e. the elliptic integrals of the first, second and third kind).

Elliptic integrals have many applications, for example in mathematics and physics

 - arc-length of plane curves (ellipse, hyperbola, Bernoulli's lemniscate)
 - surface area of ellipsoid in 3-dimensional space
 - electric and magnetic field associated with ellipsoid
 - periodicity of an anharmonic oscillator
 - mutual inductance of coaxial circles
 - age of the universe in the Friedman model
 - [Euler's three-body problem](http://en.wikipedia.org/wiki/Euler's_three-body_problem#Mathematical_solutions)
 - [and many others](http://arxiv.org/find/all/1/all:+AND+elliptic+integral/0/1/0/all/0/1).

# Elliptic Integral of the First Kind

## Definition

There are known different ways of noting the elliptic integral of the first kind, we follow here the Abramowitz et al notions.

http://latex.codecogs.com/gif.latex?\begin{aligned}F(\varphi\backslash\alpha)=F(\varphi|m)=\int^\varphi_0\frac{d\theta}{\sqrt{1-\sin^2\alpha\sin^2\theta}}=\int_0^x\frac{dt}{\sqrt{(1-t^2)(1-mt^2)}}=\int_0^udw=u\end{aligned}%.png

http://latex.codecogs.com/gif.latex?\mbox{where\,}m=\sin^2\alpha%.png


where `m=sin`^2^`alpha` and `m` is the parameter, `alpha` is the modular angle.














## Numerical evaluation

Numerical evaluation of the elliptic integrals are based on various expansion in terms of elementary and transcendental functions.

[The On-Line Encyclopedia of Integer Sequences](http://www.research.att.com/~njas/sequences/) 

### The arithmetic-geometric mean

The arithmetic-geometric mean (AGM) of two positive real numbers x and y is defined as follows:
First compute the arithmetic mean of x and y and call it a1. Next compute the geometric mean of x and y and call it g1; this is the square root of the product xy:


Then iterate this operation with a1 taking the place of x and g1 taking the place of y. In this way, two sequences (an) and (gn) are defined:


These two sequences converge to the same number, which is the arithmetic-geometric mean of x and y; it is denoted by M(x, y), or sometimes by agm(x, y).
This can be used for algorithmic purposes as in the AGM method.

### Landen's Transformation

[Landen's transformation](http://en.wikipedia.org/wiki/Landen's_transformation), independently rediscovered by Gauss, is a mapping of the parameters of an elliptic integral, which leaves the value of the integral unchanged.

### Reduction and Computation of Elliptic Integrals by B. Carlson

The conventional methods for computing elliptic integrals are Gauss and Landen transformations, which converge quadratically and work well for elliptic integrals of the first and second kinds. Unfortunately they suffer from loss of significant digits for the third kind. Carlson's algorithm, by contrast, provides a unified method for all three kinds of elliptic integrals with satisfactory precisions.


<wiki:comment>

# Elliptic Integral of the Second Kind

Elliptic integral of the second kind ... is to appear

## Definition

## Numerical evaluation

[Sequence A120362](http://www.research.att.com/~njas/sequences/A120362) represents 	 the numerators of bivariate Taylor expansion of the incomplete elliptic integral of the second kind.

# Elliptic Integral of the Third Kind

Elliptic integral of the first kind ...

## Definition

## Numerical evaluation

# Examples

## M. Abramowitz and I.A. Stegun Mathematical tables generation for elliptic integrals

## Conversion between differently defined elliptic integrals

## Elliptic integrals and derivatives

# Problem Solving with Elliptic Integrals

## Period of the physical pendulum

## Complex arguments and visualisations

</wiki:comment>

# References

  1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](http://www.math.ucla.edu/~cbm/aands/)" Dover Publications", 1965, Ch. 17.1 - 17.6.
  1. D. F. Lawden, "[Elliptic Functions and Applications](http://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989
  1. S. Zhang, J. Jin, "[Computation of Special Functions](http://jin.ece.uiuc.edu/specfunc.html)" (Wiley, 1996).
  1. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](http://nvl.nist.gov/pub/nistpubs/jres/107/5/j75car.pdf)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.

# Wiki references

  1. [Landen's transformation](http://en.wikipedia.org/wiki/Landen's_transformation): Landen's-Gauss transformation used in the elliptic integral numerical evaluations.
  1. [Elliptic integral](http://en.wikipedia.org/wiki/Elliptic_integral): the wiki community  article on the theoretical aspects of elliptic intergral
# External Links

  - [ALGLIB](http://www.alglib.net/) - Cross-platform numerical analysis and data processing library.
  - [Boost C++ Library](http://www.boost.org/doc/libs/1_36_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/ellint/ellint_intro.htmlhttp://www.boost.org/doc/libs/1_36_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/ellint/ellint_intro.html) - Elliptic Integral Overview.