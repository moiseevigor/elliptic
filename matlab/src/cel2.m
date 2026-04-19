function C = cel2(kc, a, b)
%CEL2  Bulirsch's complete elliptic integral of the second kind (generalised).
%
%   C = CEL2(KC, A, B) evaluates cel2(kc, a, b) = cel(kc, 1, a, b).
%
%   Special cases:
%       cel2(√(1-m), 1, 1-m) = E(m)    [complete 2nd kind]
%       cel2(√(1-m), 1, 0)   = B(m)    [associate B]
%       cel2(√(1-m), 0, 1)   = D(m)    [associate D]
%
%   KC, A, B may be scalars or arrays of the same size.
%
%   Reference:
%   [1] R. Bulirsch, "Numerical computation of elliptic integrals and
%       elliptic functions," Numer. Math. 7 (1965), 78–90.

if nargin < 3, error('cel2: requires three arguments (kc, a, b).'); end
C = cel(kc, 1, a, b);
end
