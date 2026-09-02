function S = weierstrassSigma(z, e1, e2, e3)
%WEIERSTRASSSIGMA  Weierstrass sigma function σ(z; e1, e2, e3).
%   S = WEIERSTRASSSIGMA(Z, E1, E2, E3) evaluates the Weierstrass sigma
%   function, which is entire, odd, and satisfies
%       σ'(z)/σ(z) = ζ(z)   (logarithmic derivative equals Weierstrass ζ)
%       σ(0) = 0,  σ'(0) = 1.
%
%   NOTE: this is the Weierstrass σ, NOT the Riemann ξ or σ functions.
%
%   Algorithm (closed theta form, DLMF 23.6.9 with 23.6.8):
%       ω1 = K(m)/√(e1-e3),  m = (e2-e3)/(e1-e3),  q = exp(-π·K(1-m)/K(m))
%       η1 = -π²/(12ω1) · θ1'''(0,q)/θ1'(0,q)
%       σ(z) = (2ω1/π) · exp(η1·z²/(2ω1)) · θ1(v,q)/θ1'(0,q),  v = πz/(2ω1)
%   The theta series converges geometrically; no quadrature is involved,
%   and every lattice zero and sign change comes out of θ1 itself.
%
%   Validity: σ(z) grows quasi-exponentially; double overflows eventually
%   through the exp(η1·z²/(2ω1)) factor.
%
%   At z = 0:  S = 0 exactly (θ1(0) = 0).
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
% Closed theta form (DLMF 23.6.9 with 23.6.8):
%
%   sigma(z) = 2*omega1/pi * exp(eta1*z^2/(2*omega1)) * theta1(v,q)/theta1'(0,q)
%   eta1     = -pi^2/(12*omega1) * theta1'''(0,q)/theta1'(0,q)
%   v = pi*z/(2*omega1),  q = exp(-pi*K(1-m)/K(m)),  m = (e2-e3)/(e1-e3)
%
% No quadrature.  The previous form integrated log(sigma) through the
% zeta pole at 2*omega1, so it was ~1e-8 inside the first period and
% catastrophically wrong (magnitude AND sign) for |z| > 2*omega1; the
% theta form is entire and carries every lattice zero and sign change.

m_param  = (e2 - e3) ./ (e1 - e3);
mp_param = (e1 - e2) ./ (e1 - e3);    % 1-m without cancellation
one = ones(size(m_param));  zed = zeros(size(m_param));
KK  = carlsonRF(zed, mp_param, one);
KKp = carlsonRF(zed, m_param,  one);
omega1 = KK ./ sqrt(e1 - e3);
q = exp(-pi .* KKp ./ KK);
v = pi .* z ./ (2 .* omega1);

[th1, ~, th1p0, th1ppp0] = weierZ_theta1_local(v, q);

eta1 = -pi^2 ./ (12 .* omega1) .* th1ppp0 ./ th1p0;
S = 2 .* omega1 ./ pi .* exp(eta1 .* z.^2 ./ (2 .* omega1)) .* th1 ./ th1p0;


% -----------------------------------------------------------------------
function [th1, th1p, th1p0, th1ppp0] = weierZ_theta1_local(v, q)
% theta1(v,q), v-derivative, theta1'(0), theta1'''(0).  The common factor
% 2 of A&S 16.27 is dropped from all four series: only the ratios
% theta1/theta1'(0) and theta1'''(0)/theta1'(0) are ever used.
qmax = max([q(:); 0]);
if qmax > 0
    nT = min(30, max(2, ceil(sqrt(abs(log(eps) ./ log(qmax))))));
else
    nT = 1;
end
th1 = zeros(size(v));  th1p = th1;
th1p0 = zeros(size(v)); th1ppp0 = th1p0;
for n = 0:nT
    qq = (-1)^n .* q.^((n+0.5)^2);
    k  = 2*n + 1;
    th1     = th1     + qq .* sin(k .* v);
    th1p    = th1p    + qq .* k .* cos(k .* v);
    th1p0   = th1p0   + qq .* k;
    th1ppp0 = th1ppp0 - qq .* k^3;
end


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
