function [a,b,c,n] = inv_agm(a0,b0,c0,tol)
% INV_AGM calculates the Inverse Artihmetic Geometric Mean of A and B (see [1]). 
% The function is used by routines ELLIPJ and ELLIPTIC12.
%
%   [A,B,C,N] = AGM(A0,B0,C0,TOL) carry out the process of the arithmetic geometric 
%   mean, starting with a given positive numbers triple (A0,B0,C0) and returns in 
%   (A,B,C) the generated sequence. N is a number of steps (returns in the type uint32).
%
%   The general scheme of the process:
%       A(i) = 1/2*( A(i-1)+B(i-1) );     A(0) = A0;
%       B(i) = sqrt( A(i-1)*B(i-1) );     B(0) = B0;
%       C(i) = 1/2*( A(i-1)+B(i-1) );     C(0) = C0;
%   Stop at the N-th step when A(N) = B(N), i.e., when C(N) = 0.
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC12, ELLIPTIC3, THETA.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 17.1 - 17.6.


if nargin<4, tol = eps; end
if nargin<3, error('Not enough input arguments.'); end

% pre-allocate space and augment if needed
chunk = 8; mmax = numel(a0);
a = zeros(chunk,mmax);
c = a;
b = a;
a(1,:) = a0;
b(1,:) = b0;
c(1,:) = c0;
n = uint32( zeros(1,mmax) );
i = 1;
while any(abs(c(i,:)) > tol)
    i = i + 1;
    if i > size(a,1)
      a = [a; zeros(chunk,mmax)];
      b = [b; zeros(chunk,mmax)];
      c = [c; zeros(chunk,mmax)];
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
[a, b, c]

% END FUNCTION AGM()