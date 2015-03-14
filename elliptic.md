

# Elliptic Functions #

**[The Jacobi's elliptic functions](http://code.google.com/p/elliptic/wiki/JacobiEllipticFunctions)** are a set of basic elliptic functions, and auxiliary theta functions, that have historical importance with also many features that show up important structure, and have direct relevance to some applications (e.g. the equation of the pendulum). They also have useful analogies to the functions of trigonometry, as indicated by the matching notation `SN` for `SIN`. They are not the simplest way to develop a general theory, as now seen: that can be said for the Weierstrass elliptic functions. They are not, however, outmoded. They were introduced by Carl Gustav Jakob Jacobi, around 1830.

**[Theta functions](http://code.google.com/p/elliptic/wiki/ThetaFunctions)** are special functions of several complex variables. They are important in several areas, including the theories of abelian varieties and moduli spaces, and of quadratic forms. They have also been applied to soliton theory. When generalized to a Grassmann algebra, they also appear in quantum field theory, specifically string theory and D-branes.

## ELLIPJ: Jacobi's elliptic functions ##

[ELLIPJ](http://code.google.com/p/elliptic/source/browse/trunk/ellipj.m) evaluates the [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi's amplitude](http://mathworld.wolfram.com/JacobiAmplitude.html).

`[Sn,Cn,Dn,Am] = ELLIPJ(U,M)` returns the values of the Jacobi elliptic functions `SN`, `CN`, `DN` and `AM` evaluated for corresponding elements of argument U and parameter M.  The arrays U and M must be of the same size (or either can be scalar).  As currently implemented, M is limited to `0 <= M <= 1`.

**General definition:**
```
u = Integral(1/sqrt(1-m^2*sin(theta)^2), 0, phi);
Sn(u) = sin(phi);
Cn(u) = cos(phi);
Dn(u) = sqrt(1-m^2*sin(phi)^2);
```

_Depends on_  `AGM`, `ELLIPKE`.<br>
<i>Used by</i> <code>THETA</code>.<br>
<i>See also</i> <code>ELLIPKE</code>.<br>
<br>
<br>
<h2>ELLIPJI: Jacobi's elliptic functions of the complex argument</h2>

<a href='http://code.google.com/p/elliptic/source/browse/trunk/ellipji.m'>ELLIPJI</a> evaluates the Jacobi elliptic functions of complex phase <code>U</code>.<br>
<br>
<code>[Sni,Cni,Dni] = ELLIPJ(U,M)</code> returns the values of the Jacobi elliptic functions <code>SNI</code>, <code>CNI</code> and <code>DNI</code> evaluated for corresponding  elements of argument <code>U</code> and parameter <code>M</code>. The arrays <code>U</code> and <code>M</code> must  be of the same size (or either can be scalar).  As currently implemented, <code>M</code> is real and limited to <code>0 &lt;= M &lt;= 1</code>.<br>
<br>
<br>
<b>Example:<br>
<pre><code>[phi1,phi2] = meshgrid(-pi:3/20:pi, -pi:3/20:pi);<br>
phi = phi1 + phi2*i;<br>
[Sni,Cni,Dni] = ellipji(phi, 0.99);<br>
</code></pre></b>

<i>Depends on</i> <code>AGM</code>, <code>ELLIPJ</code>, <code>ELLIPKE</code><br>
<i>See also</i> <code>ELLIPTIC12</code>, <code>ELLIPTIC12I</code>

<h2>JACOBITHETAETA: Jacobi's Theta and Eta Functions</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/jacobiThetaEta.m'>JACOBITHETAETA</a> evaluates Jacobi's theta and eta functions.<br>
<br>
<code>[Th, H] = JACOBITHETAETA(U,M)</code> returns the values of the Jacobi's theta and eta elliptic functions <code>TH</code> and <code>H</code> evaluated for corresponding elements of argument <code>U</code> and parameter <code>M</code>.  The arrays <code>U</code> and <code>M</code> must be the same size (or either can be scalar).  As currently implemented, <code>M</code> is real and limited to <code>0 &lt;= M &lt;= 1</code>.<br>
<br>
<b>Example:<br>
<pre><code>[phi,alpha] = meshgrid(0:5:90, 0:2:90);                  <br>
[Th, H] = jacobiThetaEta(pi/180*phi, sin(pi/180*alpha).^2);  <br>
</code></pre></b>

<i>Depends on</i> <code>AGM</code>, <code>ELLIPJ</code>, <code>ELLIPKE</code><br>
<i>See also</i> <code>ELLIPTIC12</code>, <code>ELLIPTIC12I</code>, <code>THETA</code>

<h2>THETA: Theta Functions of Four Types</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/theta.m'>THETA</a> evaluates theta functions of four types.<br>
<br>
<code>Th = THETA(TYPE,V,M)</code> returns values of theta functions<br>
evaluated for corresponding values of argument <code>V</code> and parameter <code>M</code>. <code>TYPE</code> is a type of the theta function, there are four numbered types. The arrays <code>V</code> and <code>M</code> must be the same size (or either can be scalar). As currently implemented, <code>M</code> is limited to <code>0 &lt;= M &lt;= 1</code>.<br>
<br>
<code>Th = THETA(TYPE,V,M,TOL)</code> computes the theta and eta elliptic functions to the accuracy <code>TOL</code> instead of the default <code>TOL = EPS</code>.<br>
<br>
The parameter <code>M</code> is related to the nome <code>Q</code> as <code>Q = exp(-pi*K(1-M)/K(M))</code>. Some definitions of the Jacobi's elliptic functions use the modulus <code>k</code> instead of the parameter <code>m</code>.  They are related by <code>m = k^2</code>.<br>
<br>
<b>Example:<br>
<pre><code>[phi,alpha] = meshgrid(0:5:90, 0:2:90);                  <br>
Th1 = theta(1, pi/180*phi, sin(pi/180*alpha).^2);  <br>
Th2 = theta(2, pi/180*phi, sin(pi/180*alpha).^2);  <br>
Th3 = theta(3, pi/180*phi, sin(pi/180*alpha).^2);  <br>
Th4 = theta(4, pi/180*phi, sin(pi/180*alpha).^2);  <br>
</code></pre></b>

<i>Depends on</i> <code>AGM</code>, <code>ELLIPJ</code>, <code>ELLIPKE</code>, <code>JACOBITHETAETA</code><br>
<i>See also</i> <code>ELLIPTIC12</code>, <code>ELLIPTIC12I</code>

<h1>Elliptic Integrals</h1>

<b><a href='http://code.google.com/p/elliptic/wiki/EllipticIntegrals'>Elliptic integrals</a></b> originally arose in connection with the problem of giving the arc length of an ellipse. They were first studied by Giulio Fagnano and Leonhard Euler. Modern mathematics defines an elliptic integral as any function f which can be expressed in the form<br>
<pre><code>f(x) = Integral(R(t,P(t), c, x)dt,<br>
</code></pre>
where <code>R</code> is a rational function of its two arguments, <code>P</code> is the square root of a polynomial of degree <code>3</code> or <code>4</code> with no repeated roots, and <code>c</code> is a constant.<br>
In general, elliptic integrals cannot be expressed in terms of elementary functions. Exceptions to this general rule are when <code>P</code> has repeated roots, or when <code>R(x,y)</code> contains no odd powers of <code>y</code>. However, with the appropriate reduction formula, every elliptic integral can be brought into a form that involves integrals over rational functions and the three canonical forms (i.e. the elliptic integrals of the first, second and third kind).<br>
<br>
<h2>ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/elliptic12.m'>ELLIPTIC12</a> evaluates the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function.<br>
<br>
<code>[F,E,Z] = ELLIPTIC12(U,M,TOL)</code> uses the method of the Arithmetic-Geometric Mean and Descending Landen Transformation described in <a href='http://code.google.com/p/elliptic/#References'>1</a> Ch. 17.6, to determine the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function (see <a href='http://code.google.com/p/elliptic/#References'>1</a>, <a href='http://code.google.com/p/elliptic/#References'>2</a>).<br>
<br>
<b>General definition:</b>
<pre><code>F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);<br>
E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);<br>
Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).<br>
</code></pre>


Tables generating code (see <a href='http://code.google.com/p/elliptic/wiki/elliptic#References'>1</a>, pp. 613-621):<br>
<pre><code>[phi,alpha] = meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees<br>
[F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals<br>
</code></pre>

<i>Depends on</i> <code>AGM</code><br>
<i>See also</i> <code>ELLIPKE</code>, <code>ELLIPJ</code>, <code>ELLIPTIC12I</code>, <code>ELLIPTIC3</code>, <code>THETA</code>.<br>
<br>
<h2>ELLIPTIC12I: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function of the complex argument</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/elliptic12i.m'>ELLIPTIC12i</a> evaluates the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function for the complex value of phase <code>U</code>. Parameter <code>M</code> must be in the range <code>0 &lt;= M &lt;= 1</code>.<br>
<br>
<code>[Fi,Ei,Zi] = ELLIPTIC12i(U,M,TOL)</code> where <code>U</code> is a complex phase in radians, <code>M</code> is the real parameter and <code>TOL</code> is the tolerance (optional). Default value for the tolerance is <code>eps = 2.220e-16</code>.<br>
<br>
<code>ELLIPTIC12i</code> uses the function <code>ELLIPTIC12</code> to evaluate the values of corresponding integrals.<br>
<br>
<b>Example:<br>
<pre><code>[phi1,phi2] = meshgrid(-2*pi:3/20:2*pi, -2*pi:3/20:2*pi);<br>
phi = phi1 + phi2*i;<br>
[Fi,Ei,Zi] = elliptic12i(phi, 0.5);<br>
</code></pre></b>

<i>Depends on</i> <code>ELLIPTIC12</code>, <code>AGM</code><br>
<i>See also</i> <code>ELLIPKE</code>, <code>ELLIPJ</code>, <code>ELLIPTIC3</code>, <code>THETA</code>.<br>
<br>
<h2>ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/elliptic3.m'>ELLIPTIC3</a> evaluates incomplete elliptic integral of the third kind <code>Pi = ELLIPTIC3(U,M,C)</code> where <code>U</code> is a phase in radians, <code>0 &lt; M &lt; 1</code> is the module and <code>0 &lt; C &lt; 1</code> is a parameter.<br>
<br>
<code>ELLIPTIC3</code> uses Gauss-Legendre 10 points quadrature template described in <code>[3]</code> to determine the value of the Incomplete Elliptic Integral of the Third Kind (see <code>[1, 2]</code>).<br>
<br>
<b>General definition:<br>
<pre><code>Pi(u,m,c) = int(1/((1 - c*sin(t)^2)*sqrt(1 - m*sin(t)^2)), t=0..u)<br>
</code></pre></b>

Tables generating code (<a href='http://code.google.com/p/elliptic/wiki/elliptic#References'>1</a>, pp. 625-626):<br>
<pre><code>[phi,alpha,c] = meshgrid(0:15:90, 0:15:90, 0:0.1:1);<br>
Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals<br>
</code></pre>


<h2>ELLIPTIC123: Complete and Incomplete Elliptic Integrals of the First, Second, and Third Kind</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/elliptic123.m'>ELLIPTIC123</a> is a wrapper around the different elliptic integral functions, providing a unified interface and greater range of input parameters. (Unlike ELLIPKE, ELLIPTIC12 and ELLIPTIC3, which all require a phase between zero and pi/2 and a parameter between zero and one.)<br>
<br>
<code>[F,E] = ELLIPTIC123(m)</code> — complete Elliptic Integrals of the first and second kind.<br>
<br>
<code>[F,E] = ELLIPTIC123(b,m)</code> — incomplete Elliptic Integrals of the first and second kind.<br>
<br>
<code>[F,E,PI] = ELLIPTIC123(m,n)</code> — complete Elliptic Integrals of the first to third kind.<br>
<br>
<code>[F,E,PI] = ELLIPTIC123(b,m,n)</code> — incomplete Elliptic Integrals of the first to third kind.<br>
<br>
The order of the input arguments has been chosen to be consistent with the pre-existing <code>elliptic12</code> and <code>elliptic3</code> functions.<br>
<br>
This function is still under development and its results are not always well-defined or even able to be calculated (especially for the third elliptic integral with n>1). Please see the documentation for further details.<br>
<br>
<h2>INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/inverselliptic2.m'>INVERSELLIPTIC2</a> evaluates the value of the INVERSE Incomplete Elliptic Integrals of the Second Kind.<br>
<br>
INVERSELLIPTIC2 uses the method described by Boyd J. P. to determine the value of the inverse Incomplete Elliptic Integrals of the Second Kind using the “Empirical” initialization to the Newton’s iteration method <a href='http://code.google.com/p/elliptic/wiki/elliptic#References'>7</a>.<br>
<br>
Elliptic integral of the second kind:<br>
<pre><code>E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);<br>
</code></pre>

“Empirical” initialization <a href='http://code.google.com/p/elliptic/wiki/elliptic#References'>7</a>:<br>
<pre><code>T0(z,m) = pi/2 + sqrt(r)/(theta − pi/2)<br>
</code></pre>

where<br>
<pre><code>z \in [−E(pi/2,m), E(pi/2,m)]x[0, 1] - value of the entire parameter space<br>
r = sqrt((1-m)^2 + zeta^2)<br>
zeta = 1 - z/E(pi/2,m)<br>
theta = atan((1 - m)/zeta)<br>
</code></pre>

Example:<br>
<pre><code>% modulus and phase in degrees<br>
[phi,alpha] = meshgrid(0:5:90, 0:2:90);<br>
% values of integrals<br>
[F,E] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);<br>
% values of inverse <br>
invE = inverselliptic2(E, sin(pi/180*alpha).^2);<br>
% the difference between phase phi and invE should close to zero<br>
phi - invE * 180/pi<br>
</code></pre>


<h1>Weierstrass's elliptic functions (in development)</h1>

<font color='red'>!IN DEVELOPMENT, help needed!</font> See <a href='http://groups.google.com/group/pelliptic/browse_frm/thread/f1eeff5553d3533b'>forum discussion</a>.<br>
<br>
<b><a href='http://code.google.com/p/elliptic/wiki/WeierstrassEllipticFunctions'>Weierstrass's elliptic functions</a></b> are elliptic functions that take a particularly simple form (cf Jacobi's elliptic functions); they are named for Karl Weierstrass. This class of functions are also referred to as p-functions and generally written using the symbol ℘ (a stylised letter p called Weierstrass p).<br>
<br>
The Weierstrass elliptic function can be defined in three closely related ways, each of which possesses certain advantages. One is as a function of a complex variable z and a lattice Λ in the complex plane. Another is in terms of z and two complex numbers ω1 and ω2 defining a pair of generators, or periods, for the lattice. The third is in terms z and of a modulus τ in the upper half-plane. This is related to the previous definition by τ = ω2 / ω1, which by the conventional choice on the pair of periods is in the upper half-plane. Using this approach, for fixed z the Weierstrass functions become modular functions of τ.<br>
<br>
<br>
<h1>Elliptic Related Functions</h1>

<h2>AGM: Artihmetic Geometric Mean</h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/agm.m'>AGM</a> calculates the <a href='http://en.wikipedia.org/wiki/Arithmetic-geometric_mean'>Artihmetic Geometric Mean</a> of <code>A</code> and <code>B</code> (see <a href='http://code.google.com/p/elliptic/#References'>1</a>).<br>
<br>
<code>[A,B,C,N] = AGM(A0,B0,C0,TOL)</code> carry out the process of the arithmetic geometric mean, starting with a given positive numbers triple <code>(A0, B0, C0)</code> and returns in<br>
<code>(A, B, C)</code> the generated sequence. <code>N</code> is a number of steps (returns in the value<code>uint32</code>).<br>
<br>
The general scheme of the procedure:<br>
<pre><code>A(i) = 1/2*( A(i-1)+B(i-1) );     A(0) = A0;<br>
B(i) = sqrt( A(i-1)*B(i-1) );     B(0) = B0;<br>
C(i) = 1/2*( A(i-1)+B(i-1) );     C(0) = C0;<br>
</code></pre>
Stop at the <code>N</code>-th step when <code>A(N) = B(N)</code>, i.e., when <code>C(N) = 0</code>.<br>
<br>
<i>Used by</i>  <code>ELLIPJ</code> and <code>ELLIPTIC12</code>.<br><i>See also</i> <code>ELLIPKE</code>, <code>ELLIPTIC3</code>, <code>THETA</code>.<br>
<br>
<h2>NOMEQ: The Value of Nome <code>q = q(m)</code></h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/nomeq.m'>NOMEQ</a> gives the value of Nome <code>q = q(m)</code>.<br>
<br>
Nome <code>Q = nomeq(M,TOL)</code>, where <code>0&lt;=M&lt;=1</code> is the module and <code>TOL</code> is the tolerance (optional). Default value for the tolerance is <code>eps = 2.220e-16</code>.<br>
<br>
<i>Used by</i>  <code>ELLIPJ</code>.<br>
<i>Depends on</i> <code>ELLIPKE</code><br>
<i>See also</i> <code>ELLIPTIC12I</code>, <code>ELLIPTIC3</code>, <code>THETA</code>.<br>
<br>
<h2>INVERSENOMEQ: The Value of Nome <code>m = m(q)</code></h2>
<a href='http://code.google.com/p/elliptic/source/browse/trunk/inversenomeq.m'>INVERSENOMEQ</a> gives the value of Nome <code>m = m(q)</code>.<br>
<br>
<code>M = inversenomeq(q)</code>, where <code>Q</code> is the Nome of q-series.<br>
<br>
<b>WARNING</b>. The function <code>INVERSENOMEQ</code> does not return correct values of <code>M</code> for <code>Q &gt; 0.6</code>, because of computer precision limitation. The function <code>NomeQ(m)</code> has an essential singularity at <code>M = 1</code>, so it cannot be inverted at this point and actually it is very hard to find and inverse in the neigborhood also.<br>
<br>
More preciesly:<br>
<pre><code>nomeq(1) = 1<br>
nomeq(1-eps) = 0.77548641878026<br>
</code></pre>

<b>Example:<br>
<pre><code>nomeq(inversenomeq([0.001 0.3 0.4 0.5 0.6 0.7 0.8]))<br>
</code></pre></b>

<i>Used by</i>  <code>ELLIPJ</code>.<br>
<i>Depends on</i> <code>ELLIPKE</code><br>
<i>See also</i> <code>ELLIPTIC12I</code>, <code>ELLIPTIC3</code>, <code>THETA</code>.<br>
<br>
<br>
<h1>References</h1>
<ol><li><a href='http://dlmf.nist.gov/'>NIST Digital Library of Mathematical Functions</a>
</li><li>M. Abramowitz and I.A. Stegun, "<a href='http://www.nrbook.com/abramowitz_and_stegun/'>Handbook of Mathematical Functions</a>" Dover Publications", 1965, Ch. 17.1 - 17.6.<br>
</li><li>D. F. Lawden, "<a href='http://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659'>Elliptic Functions and Applications</a>" Springer-Verlag, vol. 80, 1989<br>
</li><li>S. Zhang, J. Jin, "<a href='http://jin.ece.uiuc.edu/specfunc.html'>Computation of Special Functions</a>" (Wiley, 1996).<br>
</li><li>B. Carlson, "<a href='http://nvl.nist.gov/pub/nistpubs/jres/107/5/j75car.pdf'>Three Improvements in Reduction and Computation of Elliptic Integrals</a>", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.<br>
</li><li>N. H. Abel, "<a href='http://mathdl.maa.org/mathDL/46/?pa=content&sa=viewDocument&nodeId=1557'>Studies on Elliptic Functions</a>", english translation from french by Marcus Emmanuel Barnes. Original "Recherches sur les fonctions elliptiques", Journal fr die reine und angewandte Mathematik, Vol. 2, 1827. pp. 101-181.<br>
</li><li>B. C. Berndt, H. H. Chan, S.-S. Huang, "<a href='http://www.math.uiuc.edu/~berndt/articles/ellipticintegral.pdf'>Incomplete Elliptic Integrals in Ramanujan's Lost Notebook</a>" in q-series from a Contemporary Perspective, M. E. H. Ismail and D. Stanton, eds., Amer. Math. Soc., 2000, pp. 79-126.<br>
</li><li>J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)