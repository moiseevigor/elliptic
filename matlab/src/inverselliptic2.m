function invE = inverselliptic2(E,m,tol)
% INVERSELLIPTIC2 evaluates the value of the INVERSE Incomplete Elliptic Integrals
% of the Second Kind.
%
%   INVE = INVERSELLIPTIC2(E,M,TOL) where E is a value of the integral to
%   be inverted, 0<M<1 is the module and TOL is the tolerance (optional).
%   Default value for the tolerance is eps = 2.220e-16.
%
%   INVERSELLIPTIC2 uses the method described by Boyd J. P.
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
%
%   Example:
%       % modulus and phase in degrees
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);
%       % values of integrals
%       [F,E] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
%       % values of inverse
%       invE = inverselliptic2(E, sin(pi/180*alpha).^2);
%       % the difference between phase phi and invE should close to zero
%       phi - invE * 180/pi
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

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(E) || ~isreal(m)
    error('Input arguments must be real.');
end

if length(m)==1, m = m(ones(size(E))); end
if length(E)==1, E = E(ones(size(m))); end
if ~isequal(size(m),size(E)), error('E and M must be the same size.'); end

% check whether we've been asked to evaluate the integrals for values
% smaller than eps = 2.220446049250313e-16, if so we suppose it equal zero
m(m<eps) = 0;

invE = zeros(size(E));

% make a row vector
m = m(:);
E = E(:);

if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

% inputs
z = E; mu = 1-m;

% complete integral initialization
[~,E1] = ellipke(m,tol);

% Boyd's initialisation and the Newton iteration below only converge on
% phi in [0, pi/2].  Reduce first, using
%   E(phi + k*pi | m) = E(phi | m) + 2k*E(m)      (period)
%   E(pi - phi  | m) = 2*E(m) - E(phi | m)        (reflection)
% Oddness first: folding a tiny negative z through 2*E1 - (z + 2*E1) lost
% all its digits (rel 1e-7 at z = -1e-9*E1).
signZ = sign(z);  z = abs(z);
twoE1 = 2*E1;
k     = floor(z./twoE1);
z_red = z - k.*twoE1;                       % in [0, 2*E1)
over  = z_red > E1;
z_red(over) = twoE1(over) - z_red(over);    % now in [0, E1]

zeta = 1 - z_red./E1;
r = sqrt(zeta.*zeta+mu.*mu);
theta = atan2(mu, z_red);

% “Empirical” initialization [1]
invE(:) = pi/2 + sqrt(r).*(theta - (pi/2));
invE(:) = min(max(invE(:), 0), pi/2);

% Newton on E(phi|m) = z_red;  dE/dphi = sqrt(1 - m sin^2 phi).
% Iterate to convergence rather than a fixed count (issue #12): four steps
% are not enough near m -> 1, where the initial guess can be off by ~1.
for iter=1:100
    [~, Ecur] = elliptic12(invE(:),m,tol);
    res = Ecur - z_red;
    if all(abs(res) <= 4*eps*max(abs(z_red), realmin)), break; end   % relative
    invE(:) = invE(:) - res./max(sqrt( 1-m.*sin(invE(:)).^2 ), 1e-15);
    invE(:) = min(max(invE(:), 0), pi/2);
end

% undo the reflection, then the period strips
invE(over) = pi - invE(over);
invE(:) = signZ .* (invE(:) + k*pi);
return;

