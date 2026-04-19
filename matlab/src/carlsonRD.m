function RD = carlsonRD(x, y, z)
%CARLSONRD  Carlson's symmetric elliptic integral of the second kind R_D(x,y,z).
%
%   RD = CARLSONRD(X, Y, Z) evaluates the degenerate case R_J(x,y,z,z):
%
%       R_D(x,y,z) = (3/2) ∫_0^∞ [(t+x)(t+y)]^(-1/2) (t+z)^(-3/2) dt
%
%   for x, y >= 0 with x+y > 0, and z > 0 (DLMF 19.16.5).
%
%   Algorithm — Carlson duplication (DLMF 19.36.2):
%     Iterate  λ = √(xy) + √(yz) + √(zx),
%              S ← S + 4^{-n} / (√zₙ · (zₙ + λₙ))
%              (x,y,z) ← (x+λ, y+λ, z+λ)/4
%     until the range is small, then evaluate a 9th-order polynomial.
%
%   Connection to Legendre form:
%       D(φ|m) = (sin³φ / 3) · R_D(cos²φ, 1−m·sin²φ, 1)
%
%   X, Y, Z may be scalars or same-size arrays.
%
%   References:
%   [1] NIST DLMF §19.16, §19.36  https://dlmf.nist.gov/19
%   [2] B.C. Carlson, "Numerical Computation of Real or Complex Elliptic
%       Integrals," Numer. Algorithms 10 (1995), 13–26.

if nargin < 3, error('carlsonRD: requires three arguments (x, y, z).'); end
if ~isreal(x) || ~isreal(y) || ~isreal(z)
    error('carlsonRD: all input arguments must be real.');
end

[x, y, z] = carlsonRD_broadcast(x, y, z);
origSize = size(x);
x = x(:).';  y = y(:).';  z = z(:).';

RD = carlsonRD_core(x, y, z);
RD = reshape(RD, origSize);


% -----------------------------------------------------------------------
function RD = carlsonRD_core(x, y, z)
%CARLSONRD_CORE  Vectorised duplication algorithm (row-vector inputs).

cr    = 0.0015;   % convergence: tighter than RF because RD is 3/2 order
S     = zeros(size(x));
fac   = ones(size(x));    % 4^{-n}

for iter = 1:30
    lam  = sqrt(x.*y) + sqrt(y.*z) + sqrt(z.*x);
    sz   = sqrt(z);
    S    = S + fac ./ (sz .* (z + lam));
    fac  = fac ./ 4;
    x    = (x + lam) ./ 4;
    y    = (y + lam) ./ 4;
    z    = (z + lam) ./ 4;
    A    = (x + y + 3.*z) ./ 5;
    rng  = max([abs(x-A); abs(y-A); abs(z-A)]);
    if rng < cr * min(A)
        break;
    end
end

A  = (x + y + 3.*z) ./ 5;
X  = (A - x) ./ A;
Y  = (A - y) ./ A;
Z  = -(X + Y) ./ 3;      % (A-z)/A = (X+Y)/3  since A=(x+y+3z)/5

E2 = X.*Y - 6.*Z.^2;
E3 = (3.*X.*Y - 8.*Z.^2) .* Z;
E4 = 3.*(X.*Y - Z.^2) .* Z.^2;
E5 = X.*Y.*Z.^3;

% Polynomial from DLMF 19.36.2 Table 19.36.1
poly = 1 - 3.*E2./14 + E3./6 + 9.*E2.^2./88 - 3.*E4./22 ...
         - 9.*E2.*E3./52 + 3.*E5./26;

RD = 3.*S + fac .* A.^(-3/2) .* poly;


% -----------------------------------------------------------------------
function [x, y, z] = carlsonRD_broadcast(x, y, z)
refSz = [1 1];
for v = {x, y, z}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(x), x = x(ones(refSz)); end
if isscalar(y), y = y(ones(refSz)); end
if isscalar(z), z = z(ones(refSz)); end
if ~(isequal(size(x),size(y)) && isequal(size(y),size(z)))
    error('carlsonRD: x, y, z must be the same size or scalar.');
end
