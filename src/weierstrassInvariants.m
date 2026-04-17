function [g2, g3, Delta] = weierstrassInvariants(e1, e2, e3)
%WEIERSTRASSINVARIANTS  Lattice invariants from the three roots.
%   [G2, G3, DELTA] = WEIERSTRASSINVARIANTS(E1, E2, E3) returns the
%   Weierstrass lattice invariants G2 and G3 and the discriminant DELTA
%   for the elliptic curve  4t^3 - G2*t - G3 = 0  whose roots are E1, E2,
%   E3.
%
%   Inputs E1 > E2 > E3 must be real and satisfy E1 + E2 + E3 = 0
%   (A&S 18.1.2).  Arrays are accepted; all inputs must be the same size
%   or scalar.
%
%   Formulas (A&S 18.1.3):
%       G2    = -4*(E1*E2 + E1*E3 + E2*E3)
%       G3    =  4*E1*E2*E3
%       DELTA = G2^3 - 27*G3^2   (non-zero for a non-degenerate lattice)
%
%   Companion function: WEIERSTRASSFROМИНVARIANTS(G2, G3) -> (E1, E2, E3)
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions", Dover, 1965, §18.1.
%   [2] NIST DLMF §23.3.

if nargin < 3, error('weierstrassInvariants: requires three arguments (e1, e2, e3).'); end
if ~isreal(e1) || ~isreal(e2) || ~isreal(e3)
    error('weierstrassInvariants: e1, e2, e3 must be real.');
end

% Scalar broadcasting
if isscalar(e1) && ~isscalar(e2), e1 = e1(ones(size(e2))); end
if isscalar(e2) && ~isscalar(e1), e2 = e2(ones(size(e1))); end
if isscalar(e3) && ~isscalar(e1), e3 = e3(ones(size(e1))); end
if ~isequal(size(e1), size(e2)) || ~isequal(size(e1), size(e3))
    error('weierstrassInvariants: e1, e2, e3 must be the same size or scalar.');
end

g2    = -4 .* (e1.*e2 + e1.*e3 + e2.*e3);   % A&S 18.1.3
g3    =  4 .* e1 .* e2 .* e3;               % A&S 18.1.3
Delta = g2.^3 - 27 .* g3.^2;               % discriminant


function [e1, e2, e3] = weierstrassFromInvariants(g2, g3)
%WEIERSTRASSFROMINVARIANTS  Three roots from lattice invariants.
%   [E1, E2, E3] = WEIERSTRASSFROMINVARIANTS(G2, G3) returns the sorted
%   real roots E1 >= E2 >= E3 of  4t^3 - G2*t - G3 = 0.
%   Inputs must be scalar.  The discriminant G2^3 - 27*G3^2 must be > 0.

if nargin < 2, error('weierstrassFromInvariants: requires g2 and g3.'); end
r = sort(real(roots([4, 0, -g2, -g3])), 'descend');
e1 = r(1); e2 = r(2); e3 = r(3);
