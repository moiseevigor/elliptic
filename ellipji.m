function [sni,cni,dni] = ellipji(u,m,tol)
%ELLIPJI Jacobi elliptic functions of complex phase u.
%   [Sni,Cni,Dni] = ELLIPJ(U,M) returns the values of the Jacobi
%   elliptic functions SNI, CNI and DNI evaluated for corresponding
%   elements of argument U and parameter M.  The arrays U and M must
%   be the same size (or either can be scalar).  As currently
%   implemented, M is real and limited to 0 <= M <= 1. 
%
%   [Sni,Cni,Dni] = ELLIPJ(U,M,TOL) computes the elliptic functions to
%   the accuracy TOL instead of the default TOL = EPS.  
%
%   Some definitions of the Jacobi elliptic functions use the modulus
%   k instead of the parameter m.  They are related by m = k^2.
%
%   Example:
%   [phi1,phi2] = meshgrid(-pi:3/20:pi, -pi:3/20:pi);
%   phi = phi1 + phi2*i;
%   [Sni,Cni,Dni] = ellipji(phi, 0.99);
%
%   See also 
%       Standard: ELLIPKE, ELLIPJ, 
%       Moiseev's package: ELLIPTIC12, ELLIPTIC12I.
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
    error('The parameter m must be real.')
end

if any(m < 0) || any(m > 1) 
    error('M must be in the range 0 <= M <= 1.'); 
end

% if the input is real, evaluate the elliptic integrals with ELLIPJ
if isreal(u)
   [sni,cni,dni] = ellipj(u,m,tol);
   return;
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u))
    error('U and M must be the same size.'); 
end

% capture memory and save the structure of input arrays
sni = zeros(size(u));
cni = sni;     
dni = sni;

% make a row vector
m = m(:).'; 
u = u(:).';

% represent u in the form u = phi + i*psi
phi = real(u);
psi = imag(u);

[s,c,d] = ellipj(phi,m,tol);
[s1,c1,d1] = ellipj(psi,1-m,tol);

% function evaluations
delta = c1.^2 + m.*s.^2.*s1.^2;
sni(:) = (s.*d1 + sqrt(-1).*c.*d.*s1.*c1)./delta;
cni(:) = (c.*c1 - sqrt(-1).*s.*d.*s1.*d1)./delta;
dni(:) = (d.*c1.*d1 - sqrt(-1).*m.*s.*c.*s1)./delta;

% END FUNCTION ELLIPJI()
