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
%CEL_CORE  Vectorised evaluation (row-vector inputs).
N = numel(kc);
C = zeros(1, N);

bad = (kc < 0);
C(bad) = NaN;
pole = (p <= 0);
C(pole) = Inf;

ok = ~bad & ~pole;
if ~any(ok), return; end

m  = 1 - kc(ok).^2;
[K, ~] = ellipke(m);
[B, D, ~] = ellipticBD(m);

pp = p(ok); aa = a(ok); bb = b(ok);

% Case p ≈ 1: use a*B + b*D
% Case p ≠ 1: use Carlson/Pi formula
p1 = abs(pp - 1) < 1e-12;
pn = ~p1;

Cv = zeros(1, sum(ok));
Cv(p1) = aa(p1) .* B(p1) + bb(p1) .* D(p1);

if any(pn)
    pn_m = m(pn); pn_p = pp(pn); pn_a = aa(pn); pn_b = bb(pn); pn_K = K(pn);
    n_val = 1 - pn_p;   % n for Π(n|m) = Π(1-p|m)

    % Compute J_complete(n|m) via Carlson at φ=π/2 (s=1, c=0, d=kc):
    %   J = (1/3) * RJ(0, kc², 1, 1-n)  where 1-n = p
    kc_pn = sqrt(1 - pn_m);
    RJ_val = carlsonRJ(zeros(size(pn_m)), kc_pn.^2, ones(size(pn_m)), pn_p);
    J_n = RJ_val ./ 3;      % J_complete(n|m) = s³/3 * RJ at s=1

    Pi_n = pn_K + n_val .* J_n;    % Π(1-p|m) = K + (1-p)*J

    % Formula (DLMF §19.25, decomposition):
    %   cel = a*K + (b - a*p)*(Pi - K)/(1-p)
    Cv(pn) = pn_a .* pn_K + (pn_b - pn_a .* pn_p) .* (Pi_n - pn_K) ./ n_val;
end

C(ok) = Cv;


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
