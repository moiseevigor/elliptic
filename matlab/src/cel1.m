function C = cel1(kc)
%CEL1  Bulirsch's complete elliptic integral of the first kind.
%
%   C = CEL1(KC) evaluates cel1(kc) = cel(kc, 1, 1, 1) = K(1 - kc²).
%
%   Equivalent to the complete elliptic integral of the first kind K(m)
%   with m = 1 - kc².  The parameter kc = √(1-m) = k' (complementary
%   modulus).
%
%   KC may be a scalar or array. Requires 0 ≤ kc ≤ 1.
%
%   Reference:
%   [1] R. Bulirsch, "Numerical computation of elliptic integrals and
%       elliptic functions," Numer. Math. 7 (1965), 78–90.

if nargin < 1, error('cel1: requires one argument kc.'); end
C = cel(kc, 1, 1, 1);
end
