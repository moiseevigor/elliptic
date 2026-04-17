function S = weierstrassSigma(z, e1, e2, e3)
%WEIERSTRASSSIGMA  Weierstrass sigma function σ(z; e1, e2, e3).
%   S = WEIERSTRASSSIGMA(Z, E1, E2, E3) evaluates the Weierstrass sigma
%   function, which is entire, odd, and satisfies
%       σ'(z)/σ(z) = ζ(z)   (logarithmic derivative equals Weierstrass ζ)
%       σ(0) = 0,  σ'(0) = 1.
%
%   NOTE: this is the Weierstrass σ, NOT the Riemann ξ or σ functions.
%
%   Algorithm:
%       log σ(z) = log z + ∫_0^z [ζ(t) - 1/t] dt
%       σ(z)     = z · exp( ∫_0^z [ζ(t) - 1/t] dt )
%
%   The integral is split at  t_split = min(0.25·ω1, |z|):
%
%     ∫_0^{t_split}: Laurent series  ζ(t)-1/t = -g2·t³/60 - g3·t⁵/140 + …
%       → analytic: -g2·t_split⁴/240 - g3·t_split⁶/840
%
%     ∫_{t_split}^z: 10-pt Gauss-Legendre on ζ(t)-1/t.
%       For t ≥ t_split the GL integral for ζ(t) avoids the near-singular
%       1/t² region of ℘, so weierstrassZeta returns accurate values.
%
%   Validity: σ(z) grows quasi-exponentially; double overflows for |z|>~4ω1.
%
%   At z = 0:  S = 0 exactly.
%
%   All input conventions match WEIERSTRASSP.  Parallel and GPU modes are
%   enabled via ELLIPTIC_CONFIG.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions", Dover, 1965, §18.3, 18.5.
%   [2] NIST DLMF §23.2.

if nargin < 4, error('weierstrassSigma: requires four arguments (z, e1, e2, e3).'); end
if ~isreal(z) || ~isreal(e1) || ~isreal(e2) || ~isreal(e3)
    error('weierstrassSigma: all input arguments must be real.');
end

[z, e1, e2, e3] = weierS_broadcast(z, e1, e2, e3);
origSize = size(z);

if any(e1(:) <= e2(:)) || any(e2(:) <= e3(:))
    error('weierstrassSigma: roots must satisfy e1 > e2 > e3.');
end

S    = zeros(origSize);
N_el = numel(z);

% GPU dispatch: sigma core calls weierstrassZeta which calls weierstrassP,
% all of which auto-dispatch to GPU when has_gpu() is true.
if has_gpu()
    S = gpu_weierstrassSigma(z, e1, e2, e3, origSize);
    return;
end

% Parallel dispatch
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    S = parallel_weierstrassSigma(z, e1, e2, e3, nWorkers, minChunk, origSize);
    return;
end

% --- Serial core ---
S(:) = weierS_core(z(:).', e1(:).', e2(:).', e3(:).');


% -----------------------------------------------------------------------
function S = weierS_core(z, e1, e2, e3)
%WEIERSCORE  Vectorised serial evaluation (row-vector inputs).
%
% σ(z) = z · exp( ∫_0^z [ζ(t) - 1/t] dt )
%
% Split the integral at t_split = min(0.25·ω1, |z|):
%   [0, t_split]   → Laurent series (no cancellation)
%   [t_split, z]   → 10-pt GL on ζ(t)-1/t  (t large enough for accurate ζ)

t_gl = [ 0.9931285991850949,  0.9639719272779138, ...
         0.9122344282513259,  0.8391169718222188, ...
         0.7463319064601508,  0.6360536807265150, ...
         0.5108670019508271,  0.3737060887154195, ...
         0.2277858511416451,  0.07652652113349734 ];
w_gl = [ 0.01761400713915212, 0.04060142980038694, ...
         0.06267204833410907, 0.08327674157670475, ...
         0.1019301198172404,  0.1181945319615184,  ...
         0.1316886384491766,  0.1420961093183820,  ...
         0.1491729864726037,  0.1527533871307258   ];

N  = numel(z);
g2 = -4 .* (e1.*e2 + e1.*e3 + e2.*e3);
g3 =  4 .* e1 .* e2 .* e3;

m_param = (e2 - e3) ./ (e1 - e3);
[KK, ~] = ellipke(m_param);
omega1  = KK ./ sqrt(e1 - e3);

% Exploit odd symmetry σ(-z) = -σ(z): work with |z|, restore sign at end.
% This avoids the GL interval [t_split, z] crossing z=0 for negative z.
sgn = sign(z);
az  = abs(z);

% t_split: Laurent below, numerical ζ above.
% min(0.25*omega1, 0.99*|z|) keeps t_split inside (0, |z|).
t_split = min(0.25 .* omega1, az .* 0.99);
t_split = max(t_split, 1e-10);

% Analytic part:  ∫_0^{t_split} [ζ(t)-1/t] dt
%   ζ(t)-1/t = -g2*t³/60 - g3*t⁵/140 + O(t⁷)   (A&S 18.5.5 integrated)
%   integral  = -g2*t⁴/240 - g3*t⁶/840
int_laurent = -g2 .* t_split.^4 / 240 - g3 .* t_split.^6 / 840;

% GL part:  ∫_{t_split}^{|z|} [ζ(t)-1/t] dt
% For t >= t_split >= 0.25*omega1, the GL integral for weierstrassZeta(t)
% is free from near-singular ℘ behaviour (smallest GL node > 0.2*omega1).
a   = t_split;
b   = az;
mid = (a + b) ./ 2;
hw  = (b - a) ./ 2;
int_gl = zeros(1, N);
for k = 1:10
    tp = mid + hw .* t_gl(k);
    tm = mid - hw .* t_gl(k);
    int_gl = int_gl + w_gl(k) .* (fsigma(tp,e1,e2,e3) + fsigma(tm,e1,e2,e3));
end
int_gl = int_gl .* hw;

log_sigma = int_laurent + int_gl;
S = sgn .* az .* exp(log_sigma);

% Exact value at z = 0
S(az < eps^(1/2)) = 0;


% -----------------------------------------------------------------------
function v = fsigma(t, e1, e2, e3)
% Integrand: ζ(t) - 1/t  (called only for t >= t_split >= 0.25*omega1,
% so weierstrassZeta is accurate — no small-t cancellation here).
v = weierstrassZeta(t, e1, e2, e3) - 1 ./ t;
% Guard: at lattice poles (shouldn't occur in normal integration range)
v(~isfinite(v)) = 0;


% -----------------------------------------------------------------------
function S = gpu_weierstrassSigma(z, e1, e2, e3, origSize)
%GPU_WEIERSTRASSSIGMA  GPU path: cascade through weierstrassZeta -> GPU.
S = reshape(weierS_core(z(:).', e1(:).', e2(:).', e3(:).'), origSize);


% -----------------------------------------------------------------------
function S = parallel_weierstrassSigma(z, e1, e2, e3, nWorkers, minChunk, origSize)
%PARALLEL_WEIERSTRASSSIGMA  Split work across parfor workers.
N       = numel(z);
z_f     = z(:).'; e1_f = e1(:).'; e2_f = e2(:).'; e3_f = e3(:).';
nChunks = min(nWorkers, ceil(N / minChunk));
chunkSz = ceil(N / nChunks);
S_cells = cell(1, nChunks);
if exist('OCTAVE_VERSION', 'builtin')
    z_ch = cell(1,nChunks); e1_ch = cell(1,nChunks);
    e2_ch = cell(1,nChunks); e3_ch = cell(1,nChunks);
    for k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        z_ch{k}  = z_f(i1:i2);  e1_ch{k} = e1_f(i1:i2);
        e2_ch{k} = e2_f(i1:i2); e3_ch{k} = e3_f(i1:i2);
    end
    S_cells = parcellfun(nWorkers, @par_worker, ...
        repmat({'weierstrassSigma'}, 1, nChunks), z_ch, e1_ch, e2_ch, e3_ch, ...
        'UniformOutput', false);
else
    parfor k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        S_cells{k} = weierstrassSigma(z_f(i1:i2), e1_f(i1:i2), e2_f(i1:i2), e3_f(i1:i2));
    end
end
S = reshape([S_cells{:}], origSize);


% -----------------------------------------------------------------------
function [z, e1, e2, e3] = weierS_broadcast(z, e1, e2, e3)
refSz = [1 1];
for x = {z, e1, e2, e3}
    if numel(x{1}) > 1, refSz = size(x{1}); break; end
end
if isscalar(z),  z  = z( ones(refSz)); end
if isscalar(e1), e1 = e1(ones(refSz)); end
if isscalar(e2), e2 = e2(ones(refSz)); end
if isscalar(e3), e3 = e3(ones(refSz)); end
if ~(isequal(size(z),size(e1)) && isequal(size(e1),size(e2)) && isequal(size(e2),size(e3)))
    error('weierstrassSigma: z, e1, e2, e3 must be the same size or scalar.');
end
