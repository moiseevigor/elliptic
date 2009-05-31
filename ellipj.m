function [sn,cn,dn,am] = ellipj(u,m,tol)
%ELLIPJ Jacobi elliptic functions and Jacobi's amplitude.
%   [Sn,Cn,Dn,Am] = ELLIPJ(U,M) returns the values of the Jacobi
%   elliptic functions SN, CN, DN and AM evaluated for corresponding
%   elements of argument U and parameter M.  The arrays U and M must
%   be the same size (or either can be scalar).  As currently
%   implemented, M is limited to 0 <= M <= 1. 
%
%   [Sn,Cn,Dn,Am] = ELLIPJ(U,M,TOL) computes the elliptic functions to
%   the accuracy TOL instead of the default TOL = EPS.  
%
%   Some definitions of the Jacobi elliptic functions use the modulus
%   k instead of the parameter m.  They are related by m = k^2.
%
%   See also ELLIPKE.

%   L. Shure 1-9-88
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.14 $  $Date: 2001/04/15 12:01:40 $
%
%   Modified by Moiseev Igor, 
%   moiseev[at]sissa.it
%   34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
%   Date: 2005/10/04
%
%   The modification of orginal script fixes the problem of slow convergence
%   of the AGM algorithm for the value of parameter M=1.


%   ELLIPJ uses the method of the arithmetic-geometric mean
%   described in [1].
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 16-17.6.


if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) | ~isreal(m)
    error('Input arguments must be real.')
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

am = zeros(size(u));
cn = zeros(size(u));
sn = cn;
dn = sn;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) | any(m > 1), 
  error('M must be in the range 0 <= M <= 1.');
end

I = uint32( find(m ~= 1 & m ~= 0) );
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
	while any(abs(c(i,:)) > tol)                                    % Arithmetic-Geometric Mean of A, B and C
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
	phin(:) = (2 .^ double(n(K))).*a(i,K).*u(I);
	while i > 1
        i = i - 1;
        in = uint32( find(n(K) >= i) );
        if ~isempty(in)
          phin(in) = 0.5*(asin(c(i+1,K(in)).*sin(phin(in))./a(i+1,K(in))) + phin(in));
        end
	end
	am(I) = phin;
	sn(I) = sin(phin);
	cn(I) = cos(phin);
	dn(I) = sqrt(1 - m(I).*sin(phin).^2);
end

% Special cases: m = {0, 1} 
m0 = find(m == 0);
am(m0) = u(m0);
sn(m0) = sin(u(m0));
cn(m0) = cos(u(m0));
dn(m0) = ones(size(m0));

m1 = find(m == 1);
am(m1) = asin(tanh(u(m1)));
sn(m1) = tanh(u(m1));
cn(m1) = sech(u(m1));
dn(m1) = sech(u(m1));
