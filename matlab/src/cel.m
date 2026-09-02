function C = cel(kc, p, a, b)
%CEL  Bulirsch's generalised complete elliptic integral.
%
%   C = CEL(KC, P, A, B) evaluates Bulirsch's generalised complete
%   elliptic integral:
%
%       cel(kc,p,a,b) = ∫₀^{π/2} (a·cos²φ + b·sin²φ) /
%                        ((cos²φ + p·sin²φ)·√(cos²φ + kc²·sin²φ)) dφ
%
%   which equals ∫₀^∞ (a+b·t²)/((1+p·t²)·√((1+t²)(1+kc²·t²))) dt
%   (via the substitution t = tan φ).
%
%   Special cases:
%
%       K(m) = cel(√(1−m), 1, 1, 1)           [complete 1st kind]
%       E(m) = cel(√(1−m), 1, 1, 1−m)         [complete 2nd kind]
%       B(m) = cel(√(1−m), 1, 1, 0)           [associate B]
%       D(m) = cel(√(1−m), 1, 0, 1)           [associate D]
%       Π(n|m) = cel(√(1−m), 1−n, 1, 1)       [complete 3rd kind]
%
%   Algorithm — Carlson symmetric forms (DLMF §19.25):
%
%       m   = 1 − kc²,   kc = √(1−m)
%       K   = K(m)
%       Pi  = Π(1−p | m) = K + (1−p)·J_complete(1−p|m)
%
%       For p = 1:
%           cel = a·B(m) + b·D(m)
%
%       For p ≠ 1  (using Π = Π(1−p|m)):
%           cel = [Π·(a − 2ap + bp) + (a−b)·p·K] / (1 − p)
%
%   KC, P, A, B may be scalars or arrays of the same size (scalar inputs
%   are broadcast to the largest array).  Requires kc ≥ 0, p > 0.
%
%   Returns NaN for kc < 0; Inf for p = 0.
%
%   References:
%   [1] R. Bulirsch, "Numerical computation of elliptic integrals and
%       elliptic functions," Numer. Math. 7 (1965), 78–90.
%   [2] NIST DLMF §19.25  https://dlmf.nist.gov/19.25

% Empty input -> empty output of the same shape (elementwise semantics; the
% size checks below would otherwise reject [] against a scalar).
if nargin >= 4 && (isempty(kc) || isempty(p) || isempty(a) || isempty(b))
    sz = size(kc);
    if isempty(p), sz = size(p); end
    if isempty(a), sz = size(a); end
    if isempty(b), sz = size(b); end
    C = zeros(sz);
    return;
end

if nargin < 4, error('cel: requires four arguments (kc, p, a, b).'); end
if ~isreal(kc) || ~isreal(p) || ~isreal(a) || ~isreal(b)
    error('cel: all arguments must be real.');
end

[kc, p, a, b] = cel_broadcast(kc, p, a, b);
origSize = size(kc);
kc = kc(:).'; p = p(:).'; a = a(:).'; b = b(:).';

C = cel_core(kc, p, a, b);
C = reshape(C, origSize);


% -----------------------------------------------------------------------
function C = cel_core(kc, p, a, b)
%CEL_CORE  Bulirsch's algorithm (Numer. Math. 13 (1969) 305, "cel"), vectorised.
%   Works directly with kc, so kc ~ 1e-9 (where 1 - kc^2 rounds to 1 and the
%   previous ellipke/ellipticBD route returned Inf or garbage: cel1(1e-9) came
%   out 2e6 instead of ln(4/kc) = 22.1) and kc > 1 (m < 0) are exact, and p < 0
%   yields the Cauchy principal value.  Quadratically convergent Landen ascent;
%   the stopping test |g - k| <= g*CA leaves an error of order CA^2.
CA = 1e-9;
N  = numel(kc);
C  = nan(1, N);
k  = abs(kc);                                   % the integral depends on kc^2 only
zero_kc = (k == 0);
% kc = 0: the integrand ~ b/(p cos) at pi/2 diverges unless b = 0; for b = 0
% the limit kc -> 0 is finite and the ascent evaluates it from realmin.
k(zero_kc & (b == 0)) = realmin;
run = ~zero_kc | (b == 0);
if any(zero_kc & (b ~= 0))
    C(zero_kc & (b ~= 0)) = sign(b(zero_kc & (b ~= 0)) ./ p(zero_kc & (b ~= 0))) .* Inf;
end
if ~any(run), return; end
k = k(run); p = p(run); a = a(run); b = b(run);
n  = numel(k);
e  = k;  em = ones(1, n);
pos = p > 0;
% p > 0
pp = sqrt(p(pos));  p(pos) = pp;  b(pos) = b(pos) ./ pp;
% p <= 0: transform to the p > 0 case (principal value); p = 0 gives Inf below
neg = ~pos;
if any(neg)
    f = k(neg).^2;  q = 1 - f;  g = 1 - p(neg);  f = f - p(neg);
    q = q .* (b(neg) - a(neg) .* p(neg));
    pn = sqrt(f ./ g);
    an = (a(neg) - b(neg)) ./ g;
    b(neg) = -q ./ (g.^2 .* pn) + an .* pn;
    a(neg) = an;  p(neg) = pn;
end
active = true(1, n);
for it = 1:60
    f = a(active);
    a(active) = a(active) + b(active) ./ p(active);
    g = e(active) ./ p(active);
    b(active) = 2 .* (b(active) + f .* g);
    p(active) = p(active) + g;
    g = em(active);
    em(active) = em(active) + k(active);
    conv = abs(g - k(active)) <= g .* CA;
    idx = find(active);
    kk = 2 .* sqrt(e(idx(~conv)));
    k(idx(~conv)) = kk;
    e(idx(~conv)) = kk .* em(idx(~conv));
    active(idx(conv)) = false;
    if ~any(active), break; end
end
C(run) = pi / 2 .* (b + a .* em) ./ (em .* (em + p));

% -----------------------------------------------------------------------
function [kc, p, a, b] = cel_broadcast(kc, p, a, b)
refSz = [1 1];
for v = {kc, p, a, b}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(kc), kc = kc(ones(refSz)); end
if isscalar(p),  p  = p(ones(refSz));  end
if isscalar(a),  a  = a(ones(refSz));  end
if isscalar(b),  b  = b(ones(refSz));  end
if ~(isequal(size(kc),size(p)) && isequal(size(p),size(a)) && isequal(size(a),size(b)))
    error('cel: kc, p, a, b must be the same size or scalar.');
end
