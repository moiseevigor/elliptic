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
%     moiseev.igor[at]gmail.com, moiseev[at]sissa.it
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
lambda = []; mu = []; 
I = [];      J  = [];

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

% finding the roots of the equation
% X^2 - (cot(phi)^2+m*sinh(psi)^2*csc(phi)^2-1+m)X - (1-m)*cot(phi)^2 = 0
b = -(cot(phi).^2 + m.*sinh(psi).^2.*csc(phi).^2-1+m);
c = -(1-m).*cot(phi).^2;

X1 = -b/2 + sqrt(b.^2/4-c);
I = find(X1>=0);

if length(I) ~= length(u)
    X2 = -b/2 - sqrt(b.^2/4-c);
    J = find(X2>=0);
end

if( ~isempty(I) ) 
    lambda(I) = acot( sqrt(X1(I)) ); 
    mu(I)     = atan( sqrt(1./m(I).*(tan(phi(I)).^2.*cot(lambda(I)).^2 - 1)) );
end
if( ~isempty(J) ) 
    lambda(J) = acot( sqrt(X2(J)) ); 
    mu(J)     = atan( sqrt(1./m(J).*(tan(phi(J)).^2.*cot(lambda(J)).^2 - 1)) );
end

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

% END FUNCTION ELLIPTIC12i()
