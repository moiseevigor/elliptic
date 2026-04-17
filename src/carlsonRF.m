function RF = carlsonRF(x, y, z)
%CARLSONRF  Carlson's symmetric elliptic integral of the first kind R_F(x,y,z).
%
%   RF = CARLSONRF(X, Y, Z) evaluates
%
%       R_F(x,y,z) = (1/2) ∫_0^∞ [(t+x)(t+y)(t+z)]^(-1/2) dt
%
%   for x, y, z >= 0 with at most one zero argument (DLMF 19.16.1).
%
%   Algorithm — Carlson duplication (DLMF 19.36.1):
%     Iterate  λ = √(xy) + √(yz) + √(zx),
%              (x,y,z) ← (x+λ, y+λ, z+λ)/4
%     until the range max|arg − mean| < cr ≈ 0.0027.
%     Then evaluate the 7th-degree polynomial in the elementary symmetric
%     functions E₂, E₃.
%
%   Connection to Legendre form:
%       F(φ|m) = sin(φ) · R_F(cos²φ, 1−m·sin²φ, 1)
%
%   X, Y, Z may be any combination of scalars or same-size arrays.
%
%   No separate parallel/GPU dispatch at this level; higher-level routines
%   own the dispatch logic.
%
%   References:
%   [1] NIST DLMF §19.16, §19.36  https://dlmf.nist.gov/19
%   [2] B.C. Carlson, "Numerical Computation of Real or Complex Elliptic
%       Integrals," Numer. Algorithms 10 (1995), 13–26.

if nargin < 3, error('carlsonRF: requires three arguments (x, y, z).'); end
if ~isreal(x) || ~isreal(y) || ~isreal(z)
    error('carlsonRF: all input arguments must be real.');
end

[x, y, z] = carlsonRF_broadcast(x, y, z);
origSize = size(x);
x = x(:).';  y = y(:).';  z = z(:).';

RF = carlsonRF_core(x, y, z);
RF = reshape(RF, origSize);


% -----------------------------------------------------------------------
function RF = carlsonRF_core(x, y, z)
%CARLSONRF_CORE  Vectorised duplication algorithm (row-vector inputs).

% Duplication convergence criterion for double precision (DLMF 19.36.1)
cr = 0.0027;

x0 = x;  y0 = y;  z0 = z;

for iter = 1:20
    lam = sqrt(x.*y) + sqrt(y.*z) + sqrt(z.*x);
    x   = (x + lam) ./ 4;
    y   = (y + lam) ./ 4;
    z   = (z + lam) ./ 4;
    A   = (x + y + z) ./ 3;
    if max(max(abs(x - A)), max(max(abs(y - A)), abs(z - A))) < cr * min(A)
        break;
    end
end

A  = (x + y + z) ./ 3;
X  = (A - x) ./ A;
Y  = (A - y) ./ A;
Z  = -X - Y;          % X+Y+Z = 0 by construction

E2 = X.*Y - Z.^2;
E3 = X.*Y.*Z;

% Polynomial fit accurate to ~1e-16 (DLMF 19.36.1 Table)
RF = A.^(-0.5) .* (1 - E2./10 + E3./14 + E2.^2./24 - 3.*E2.*E3./44);


% -----------------------------------------------------------------------
function [x, y, z] = carlsonRF_broadcast(x, y, z)
refSz = [1 1];
for v = {x, y, z}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(x), x = x(ones(refSz)); end
if isscalar(y), y = y(ones(refSz)); end
if isscalar(z), z = z(ones(refSz)); end
if ~(isequal(size(x),size(y)) && isequal(size(y),size(z)))
    error('carlsonRF: x, y, z must be the same size or scalar.');
end
