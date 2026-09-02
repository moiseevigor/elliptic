function [Fi,Ei,Zi] = elliptic12i(u,m,tol)
% ELLIPTIC12i evaluates the Incomplete Elliptic Integrals
% of the First, Second Kind and Jacobi's Zeta Function for the complex
% value of phase U. Parameter M must be in the range 0 <= M <= 1.
%
%   [Fi,Ei,Zi] = ELLIPTIC12i(U,M,TOL) where U is a complex phase in
%   radians, M is the real parameter and TOL is the tolerance (optional).
%   Default value for the tolerance is eps = 2.220e-16.
%
%   ELLIPTIC12i uses the function ELLIPTIC12 to evaluate the values of
%   corresponding integrals.
%
%   Branch convention: values follow the A&S 17.4.11 real decomposition
%   F(phi+i*psi|m) = F(lambda|m) + i*F(mu|1-m).  On the line phi = pi/2
%   (above the branch point pi/2 + i*acosh(1/sqrt(m))) this fixes
%   Re F = K(m), which diverges as m -> 1.  Other systems (Mathematica,
%   mpmath) may return values on a different sheet there; both satisfy
%   the defining differential relation.
%
%   Example:
%   [phi1,phi2] = meshgrid(-2*pi:3/20:2*pi, -2*pi:3/20:2*pi);
%   phi = phi1 + phi2*i;
%   [Fi,Ei,Zi] = elliptic12i(phi, 0.5);
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC12.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html
% Everyone is permitted to copy and distribute verbatim copies of this
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE.
%
% Copyright (C) 2007 by Moiseev Igor. All rights reserved.
% 34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
% For support, please reply to
%     moiseev.igor[at]gmail.com
%     Moiseev Igor,
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(m)
    error('The parameter M must be real.')
end

if any(m < 0) || any(m > 1)
    error('M must be in the range 0 <= M <= 1.');
end

% if the input is real, evaluate the elliptic integrals with ELLIPTIC12
% if isreal(u)
%    [Fi,Ei,Zi] = elliptic12(u,m,tol);
%    return;
% end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u))
    error('U and M must be the same size.');
end

% capture memory and save the structure of input arrays
F1 = zeros(size(u)); F2 = zeros(size(u));
E1 = F1;     E2 = F1;
Z1 = F1;     Z2 = F1;
Fi = F1;     Ei = F1;
Zi = F1;

% make a row vector
m = m(:).';
u = u(:).';

% represent u in the form u = phi + i*psi
phi = real(u);
psi = imag(u);

% to avoid singularity of COT(phi) at zero add EPS
I = find (abs(phi) < eps);
phi(I) = eps;
I = [];

% finding the roots of the equation, with X = cot(lambda)^2
% X^2 - (cot(phi)^2+m*sinh(psi)^2*csc(phi)^2-1+m)X - (1-m)*cot(phi)^2 = 0
cot2 = cot(phi).^2;
b = -(cot2 + m.*sinh(psi).^2.*csc(phi).^2-1+m);
c = -(1-m).*cot2;

% Positive root X1 = cot(lambda)^2 of X^2 + bX + c = 0 and tan(mu)^2, both
% without cancellation.  Writing X1 = cot(phi)^2 + Y, Y solves
%     Y^2 + B'Y - C' = 0,  B' = cot2 + (1-m) - m sinh^2 csc^2,
%                          C' = cot2 * m sinh^2 csc^2 >= 0,
% and A&S 17.4.11's tan(mu)^2 = (tan(phi)^2 cot(lambda)^2 - 1)/m = Y/(m cot2)
% collapses to
%     tan(mu)^2 = 2 sinh^2 csc^2 / (B' + sqrt(B'^2 + 4C'))        (B' >= 0)
%               = (|B'| + sqrt(B'^2 + 4C')) / (2 m cot2)          (B' <  0)
% -- m cancels analytically in the first form, so m -> 0 (and m = 0 exactly)
% is handled to full precision.  The old (ratio-1)/m lost sqrt(eps/m) digits
% and returned Im F = 0 for psi = 1e-9.
s2c2 = sinh(psi).^2.*csc(phi).^2;
Bp   = cot2 + (1-m) - m.*s2c2;
Cp   = cot2.*m.*s2c2;
root = sqrt(Bp.^2 + 4*Cp);
pos  = Bp >= 0;
Y      = zeros(size(Bp));  tan2mu = Y;
Y(pos)  = 2*Cp(pos)./(Bp(pos) + root(pos));
Y(~pos) = 0.5*(-Bp(~pos) + root(~pos));
tan2mu(pos)  = 2*s2c2(pos)./(Bp(pos) + root(pos));
tan2mu(~pos) = 0.5*(-Bp(~pos) + root(~pos))./(m(~pos).*cot2(~pos));
X1 = cot2 + Y;

lambda = acot( sqrt(X1) );
mu     = atan( sqrt(tan2mu) );

% change of variables taking into account periodicity ceil to the right
lambda = (-1).^floor(phi/pi*2).*lambda + pi*ceil(phi/pi-0.5+eps);
mu     = sign(psi).*real(mu);

[F1(:),E1(:)] = elliptic12(lambda, m, tol);
[F2(:),E2(:)] = elliptic12(mu, 1-m, tol);

% complex values of elliptic integral of the first kind
Fi = F1 + sqrt(-1)*F2;

% some calucation optimiziation
sin_lam = sin(lambda); cos_lam = cos(lambda);
sin_mu = sin(mu); cos_mu = cos(mu);

b1 = m.*sin_lam.*cos_lam.*sin_mu.^2.*sqrt(1-m.*sin_lam.^2);
b2 = sin_mu.*cos_mu.*(1-m.*sin_lam.^2).*sqrt(1-(1-m).*sin_mu.^2);
b3 = cos_mu.^2 + m.*sin_lam.^2.*sin_mu.^2;

% complex values of elliptic integral of the second kind
Ei(:) = (b1 + sqrt(-1)*b2)./b3;
Ei(:) = Ei(:) + E1(:) + sqrt(-1)*(-E2(:) + F2(:));

[K,Ee] = ellipke(m);
% complex values of zeta function
Zi(:) = Ei(:) - Ee(:)./K(:).*Fi(:);

% Small-m Maclaurin series (through m^2).  The A&S 17.4.11 decomposition
% loses ~sqrt(eps/m) digits as m -> 0 (0.2 absolute at m = 1e-16); the
% series is exact there and covers m = 0 itself:
%   F = u + m(u/4 - sin2u/8) + m^2(9u/64 - 3sin2u/32 + 3sin4u/256) + O(m^3)
%   E = u - m(u/4 - sin2u/8) - m^2(3u/64 -  sin2u/32 +  sin4u/256) + O(m^3)
% Valid while |m sin^2 u| is small: switch on m*max(1, e^(2|psi|)) < 1e-4,
% where the crossover error is ~2e-12 (measured against 40-digit mpmath).
m_eff = m .* max(1, exp(2*abs(psi)));
sm = find(m_eff < 1e-4);
if ~isempty(sm)
    uu = u(sm);  mm = m(sm);    % original u: phi carries the +eps nudge
    s2 = sin(2*uu);  s4 = sin(4*uu);
    Fs = uu + mm.*(uu/4 - s2/8) + mm.^2.*(9*uu/64 - 3*s2/32 + 3*s4/256);
    Es = uu - mm.*(uu/4 - s2/8) - mm.^2.*(3*uu/64 - s2/32 + s4/256);
    Fi(sm) = Fs;
    Ei(sm) = Es;
    Zi(sm) = Es - Ee(sm)./K(sm).*Fs;
end

% END FUNCTION ELLIPTIC12i()
