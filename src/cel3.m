function C = cel3(kc, p)
%CEL3  Bulirsch's complete elliptic integral of the third kind (generalised).
%
%   C = CEL3(KC, P) evaluates cel3(kc, p) = cel(kc, p, 1, 1).
%
%   Connection to Legendre's complete third kind:
%       cel3(√(1-m), 1-n) = Π(n|m)    [complete 3rd kind, sign convention]
%
%   KC, P may be scalars or arrays of the same size.
%
%   Reference:
%   [1] R. Bulirsch, "Numerical computation of elliptic integrals and
%       elliptic functions," Numer. Math. 7 (1965), 78–90.

if nargin < 2, error('cel3: requires two arguments (kc, p).'); end
C = cel(kc, p, 1, 1);
end
