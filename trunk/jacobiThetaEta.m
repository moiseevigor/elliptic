function [Th,H] = jacobiThetaEta(u,m,tol)
%JACOBITHETAETA evaluates Jacobi's theta and eta functions.
%   [Th, H] = JACOBITHETAETA(U,M) returns the values of the Jacobi's
%   theta and eta elliptic functions TH and H evaluated for corresponding
%   elements of argument U and parameter M.  The arrays U and M must
%   be the same size (or either can be scalar).  As currently
%   implemented, M is limited to 0 <= M <= 1. 
%
%   [Th, H] = JACOBITHETAETA(U,M,TOL) computes the theta and eta 
%   elliptic functions to the accuracy TOL instead of the default TOL = EPS.  
%
%   Some definitions of the Jacobi elliptic functions use the modulus
%   k instead of the parameter m.  They are related by m = k^2.
%
%   Example:
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);                  
%       [Th, H] = jacobiThetaEta(pi/180*phi, sin(pi/180*alpha).^2);  
%
%   See also 
%       Standard: ELLIPKE, ELLIPJ, 
%       Moiseev's package: ELLIPTIC12, ELLIPTIC12I, THETA.
%
%   ELLIPJ uses the method of the arithmetic-geometric mean
%   described in [1].
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 16-17.6.

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

if ~isreal(u) | ~isreal(m)
    error('Input arguments must be real.')
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

Th = zeros(size(u));
H = Th;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) | any(m > 1), 
  error('M must be in the range 0 <= M <= 1.');
end

KK = ellipke(m);
period_condition = u./KK/2-floor(u./KK/2);

I_odd = uint32( find( abs(m-1) > 10*eps & abs(m) > 10*eps & abs(period_condition - 0.5) < 10*eps ) );
% here we cheat, add some disturbance to avoid the AGM algorithm divergence
% when the inputs are the ratios of complete elliptic integrals
if ( ~isempty(I_odd) )
    u(I_odd) = u(I_odd) + 100000*eps;
    m(I_odd) = m(I_odd) + 10000*eps;
end

I = uint32( find( abs(m-1) > 10*eps & ...
                  abs(m) > 10*eps ... 
                 ) ...
           );
%                   abs(period_condition - 0.5) > 10*eps ...                % odd period
%                   abs(period_condition) > 10*eps ...                      % even period

if ~isempty(I)
    [mu,J,K] = unique(m(I));   % extracts unique values from m
    K = uint32(K);
    mumax = length(mu);

    % pre-allocate space and augment if needed
	chunk = 7;
	a = zeros(chunk,mumax);
	c = a; 
	b = a;
	a(1,:) = ones(1,mumax);
	c(1,:) = sqrt(mu);
	b(1,:) = sqrt(1-mu);
	n = uint32( zeros(1,mumax) );
	i = 1;
    
    % Arithmetic-Geometric Mean of A, B and C
    while any(abs(c(i,:)) > tol)                                    
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,mumax)];
          b = [b; zeros(2,mumax)];
          c = [c; zeros(2,mumax)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        in = uint32( find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol)) );
        if ~isempty(in)
          [mi,ni] = size(in);
          n(in) = ones(mi,ni)*(i-1);
        end
	end

    mmax = length(I);
	phin = zeros(1,mmax);
    prodth = ones(i,mmax);

    phin(:) = (2 .^ double(n(K))).*a(i,K).*u(I);
    phin_pred = phin;
	while i > 1
        i = i - 1;
        in = uint32( find(n(K) >= i) );
        if ~isempty(in)
          phin(in) = 0.5*(asin(c(i+1,K(in)).*sin(phin(in))./a(i+1,K(in))) + phin(in));
          prodth(i,in) = ( sec(2*phin(in)-phin_pred(in)) ).^(1/2^(i+1));
          if (i > 1)
              phin_pred = phin;
          end
        end
    end
    
    th_save = sqrt(2*sqrt(1-m(I)).* KK(I)/pi.* cos(phin_pred - phin)./cos(phin) ).* prod(prodth,1);
    Th(I) = th_save;
    H(I) = sqrt(sqrt(m(I))).* sin(phin).* th_save; 
end

% special values of u = (2n+1)*KK, odd periods
% I_odd = uint32( find( abs(m-1) > 10*eps & abs(m) > 10*eps & abs(period_condition - 0.5) < 10*eps ) );
% if ( ~isempty(I_odd) )
%     Th(I_odd) = 1+2*(1./(1-exp(-pi*ellipke(1-m(I_odd))./KK(I_odd)))-1);
%     Th(I_odd) = sqrt(KK(I_odd)./ellipke(1-m(I_odd)));
%     H(I_odd)  = (-1).^(floor(u(I_odd)./KK(I_odd)/2)).* sqrt(sqrt(m(I_odd))).* Th(I_odd); 
% end
    
% Special cases: m = {0, 1} 
m0 = find(abs(m) < 10*eps);

if ( ~isempty(m0) )
    Th(m0) = 1;
    H(m0)  = sqrt(sqrt(m(m0))).* sin(u(m0));
end

m1 = find(abs(m-1) < 10*eps);
if ( ~isempty(m0) )
    Th(m1) = NaN;
    H(m1)  = NaN;
end