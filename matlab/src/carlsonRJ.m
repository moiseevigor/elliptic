function RJ = carlsonRJ(x, y, z, p)
%CARLSONRJ  Carlson's symmetric elliptic integral of the third kind R_J(x,y,z,p).
%
%   RJ = CARLSONRJ(X, Y, Z, P) evaluates
%
%       R_J(x,y,z,p) = (3/2) ∫_0^∞ [(t+x)(t+y)(t+z)]^(-1/2) (t+p)^(-1) dt
%
%   for x, y, z >= 0 with x+y+z > 0, and p ≠ 0 (DLMF 19.16.2).
%
%   Algorithm — Carlson duplication (DLMF 19.36.3):
%     At each step compute
%        δ  = (p + λ)(p₀ − p)   [accumulating correction]
%        S  += 4^{-n} · R_C(1, 1 + δ/(4ⁿ p²))
%     then iterate (x,y,z,p) ← (x+λ,y+λ,z+λ,p+λ)/4.
%     Final evaluation uses a 11th-order polynomial.
%
%   Connection to Legendre form:
%       Π(n,φ|m) = sin(φ) · [R_F(c²,d²,1) + (n/3)·sin²φ · R_J(c²,d²,1,1−n·s²)]
%       where  s=sinφ, c=cosφ, d=√(1−m·s²)
%
%   X, Y, Z, P may be scalars or same-size arrays.
%
%   References:
%   [1] NIST DLMF §19.16, §19.36  https://dlmf.nist.gov/19
%   [2] B.C. Carlson, "Numerical Computation of Real or Complex Elliptic
%       Integrals," Numer. Algorithms 10 (1995), 13–26.

if nargin < 4, error('carlsonRJ: requires four arguments (x, y, z, p).'); end
if ~isreal(x) || ~isreal(y) || ~isreal(z) || ~isreal(p)
    error('carlsonRJ: all input arguments must be real.');
end

[x, y, z, p] = carlsonRJ_broadcast(x, y, z, p);
origSize = size(x);
x = x(:).';  y = y(:).';  z = z(:).';  p = p(:).';

RJ = carlsonRJ_core(x, y, z, p);
RJ((x == 0) + (y == 0) + (z == 0) >= 2) = Inf;   % diverges (DLMF 19.16.2)
RJ = reshape(RJ, origSize);


% -----------------------------------------------------------------------
function RJ = carlsonRJ_core(x, y, z, p)
%CARLSONRJ_CORE  Vectorised duplication algorithm (row-vector inputs).

cr   = 0.0015;
S    = zeros(size(x));
fac  = ones(size(x));

p0   = p;    % save original p for δ computation

% Each duplication divides the argument-ratio exponent (base 4) by one; the
% adaptive break below decides, the cap only guards pathological input.
% A cap of 30 covered ratios to ~1e16 only (RJ(1e-20,2e-20,3e-20,.5) 11% off).
% Per-element convergence: each element stops when IT has converged, so a
% value never depends on what else is in the batch (a vector-wide break
% made chunked and serial evaluations differ by an ulp).
active = true(size(x));
for iter = 1:200
    lam  = sqrt(x.*y) + sqrt(y.*z) + sqrt(z.*x);
    % R_C argument for sum term (DLMF 19.36.3)
    alpha = (p .* (sqrt(x) + sqrt(y) + sqrt(z)) + sqrt(x.*y.*z)).^2;
    beta  = p .* (p + lam).^2;
    S(active)   = S(active) + fac(active) .* carlsonRC_core(alpha(active), beta(active));
    fac(active) = fac(active) ./ 4;
    x(active) = (x(active) + lam(active)) ./ 4;
    y(active) = (y(active) + lam(active)) ./ 4;
    z(active) = (z(active) + lam(active)) ./ 4;
    p(active) = (p(active) + lam(active)) ./ 4;
    A     = (x + y + z + 2.*p) ./ 5;
    active = active & (max([abs(x-A); abs(y-A); abs(z-A); abs(p-A)], [], 1) >= cr * A);
    if ~any(active), break; end
end

A  = (x + y + z + 2.*p) ./ 5;
X  = (A - x) ./ A;
Y  = (A - y) ./ A;
Z  = (A - z) ./ A;
P  = -(X + Y + Z) ./ 2;    % (A-p)/A

E2 = X.*Y + X.*Z + Y.*Z - 3.*P.^2;
E3 = X.*Y.*Z + 2.*E2.*P + 4.*P.^3;   % DLMF 19.36.2 (was 3P^3: O(eps^4) truncation, 1e-13 relative)
E4 = (2.*X.*Y.*Z + E2.*P + 3.*P.^3) .* P;
E5 = X.*Y.*Z.*P.^2;

% Polynomial from DLMF 19.36.3 Table 19.36.1
poly = 1 - 3.*E2./14 + E3./6 + 9.*E2.^2./88 - 3.*E4./22 ...
         - 9.*E2.*E3./52 + 3.*E5./26;

RJ = 3.*S + fac .* A.^(-3/2) .* poly;


% -----------------------------------------------------------------------
function RC = carlsonRC_core(x, y)
%CARLSONRC_CORE  Shared RC implementation (duplicate to avoid circular dependency).
N  = numel(x);
RC = zeros(1, N);

pole   = (y == 0);
zero_x = (x == 0) & ~pole;
eq     = (x == y) & ~pole & ~zero_x;
gt     = (y > x) & ~pole & ~zero_x & ~eq;
lt     = (y < x) & ~pole & ~zero_x & ~eq;

RC(pole)   = Inf;
RC(zero_x) = pi ./ (2 .* sqrt(y(zero_x)));
RC(eq)     = 1 ./ sqrt(x(eq));

if any(gt)
    d      = sqrt((y(gt) - x(gt)) ./ x(gt));
    RC(gt) = atan(d) ./ sqrt(y(gt) - x(gt));
end
if any(lt)
    % log((sqrt(x)+sqrt(x-y))/sqrt(y))/sqrt(x-y) == atanh(sqrt(1-y/x))/sqrt(x-y)
    % without the 1 - sqrt(1-eps) cancellation (RC(3,1e-10) lost 8 digits).
    % ... as log1p: log((sx+sxy)/sy) = log1p(((x-y)/(sx+sy) + sxy)/sy); the
    % plain log lost 9 digits again for tiny x - y (RC(1+1e-13, 1)).
    xl = x(lt);  yl = y(lt);  sx = sqrt(xl);  sy = sqrt(yl);  sxy = sqrt(xl - yl);
    RC(lt) = log1p(((xl - yl)./(sx + sy) + sxy) ./ sy) ./ sxy;
end


% -----------------------------------------------------------------------
function [x, y, z, p] = carlsonRJ_broadcast(x, y, z, p)
refSz = [1 1];
for v = {x, y, z, p}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(x), x = x(ones(refSz)); end
if isscalar(y), y = y(ones(refSz)); end
if isscalar(z), z = z(ones(refSz)); end
if isscalar(p), p = p(ones(refSz)); end
if ~(isequal(size(x),size(y)) && isequal(size(y),size(z)) && isequal(size(z),size(p)))
    error('carlsonRJ: x, y, z, p must be the same size or scalar.');
end
