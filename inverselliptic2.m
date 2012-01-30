function invE = inverselliptic2(E,m,tol)
% INVERSELLIPTIC evaluates the value of the INVERSE Incomplete Elliptic Integrals 
% of the Second Kind.
%
%   INVE = INVERSELLIPTIC(E,M,TOL) where E is a value of the integral to 
%   be converted, 0<M<1 is the module and TOL is the tolerance (optional). 
%   Default value for the tolerance is eps = 2.220e-16.
%
%   INVERSELLIPTIC uses the method described by Boyd J. P. 
%   to determine the value of the inverse Incomplete Elliptic Integrals 
%   of the Second Kind using the “Empirical” initialization to 
%   the Newton’s iteration method [1].
%
%   NOTICE. Please pay attention to the definition of the elliptic functions
%   which follows the Abramovitz et al [2], for more theory on elliptic
%   functions please consult the Lawden book [3].
%
%   Elliptic integral of the second kind:
%
%       E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
%
%   “Empirical” initialization [1]:
%
%       T0(z,m) = pi/2 + sqrt(r)/(theta − pi/2)
%
%   where 
%       z \in [−E(pi/2,m), E(pi/2,m)]x[0, 1], value of the entire parameter space
%       r = sqrt((1-m)^2 + zeta^2)
%       zeta = 1 - z/E(pi/2,m)
%       theta = atan((1 - m)/zeta)
%
%   See also ELLIPKE, ELLIPTIC12.
%
%   References:
%   [1] J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion 
%       of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)
%   [2] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
%   [3] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989
%
% Copyright (C) 2011 by Elliptic Project. All rights reserved.

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html 
% Everyone is permitted to copy and distribute verbatim copies of this 
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE. 
%  
%   For support, please reply to 
%       moiseev.igor[at]gmail.com
%       Moiseev Igor
%
%   ELLIPTIC PROJECT: http://elliptic.googlecode.com
%   Group: 

z=E; mu=1-m;
[K,E1] = ellipke(m); 
zeta = 1 - z./E1;
r = sqrt(zeta.*zeta+mu.*mu); 
theta = atan(mu./(z+1.E-12));
invE = pi/2 + sqrt(r).*(theta - (pi/2)); 

for iter=1:4
    [F, E] = elliptic12(invE,m);
    invE = invE-(E - z)./sqrt(1-m.*sin(invE).^2 ); 
end
return;

