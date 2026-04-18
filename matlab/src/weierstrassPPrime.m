function dP = weierstrassPPrime(z, e1, e2, e3)
%WEIERSTRASSP PRIME  Derivative ℘'(z; e1, e2, e3) of the Weierstrass P-function.
%   dP = WEIERSTRASSPPRIME(Z, E1, E2, E3) evaluates the derivative of the
%   Weierstrass P-function at argument Z.
%
%   Formula (chain rule on A&S 18.9.1):
%       m   = (e2 - e3) / (e1 - e3)
%       w   = z * sqrt(e1 - e3)
%       dP  = -2 * (e1 - e3)^(3/2) * cn(w,m) * dn(w,m) / sn^3(w,m)
%
%   This is the explicit form of A&S 18.9.8.  The sign is determined
%   automatically; no sqrt-with-sign ambiguity.
%
%   Verification identity (A&S 18.1.5):
%       (℘')^2 = 4*℘^3 - g2*℘ - g3
%   where g2 = -4*(e1*e2+e1*e3+e2*e3), g3 = 4*e1*e2*e3.
%
%   At half-period ω1 = K(m)/√(e1-e3): ℘'(ω1) = 0 exactly.
%   At poles (z = 0 or lattice points): dP = Inf.
%
%   All input conventions match WEIERSTRASSP.  Parallel and GPU modes are
%   enabled via ELLIPTIC_CONFIG.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions", Dover, 1965, §18.9.

if nargin < 4, error('weierstrassPPrime: requires four arguments (z, e1, e2, e3).'); end
if ~isreal(z) || ~isreal(e1) || ~isreal(e2) || ~isreal(e3)
    error('weierstrassPPrime: all input arguments must be real.');
end

[z, e1, e2, e3] = weierPP_broadcast(z, e1, e2, e3);
origSize = size(z);

if any(e1(:) <= e2(:)) || any(e2(:) <= e3(:))
    error('weierstrassPPrime: roots must satisfy e1 > e2 > e3.');
end

dP   = zeros(origSize);
N_el = numel(z);

% GPU dispatch
if has_gpu()
    dP = gpu_weierstrassPPrime(z, e1, e2, e3, origSize);
    return;
end

% Parallel dispatch
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    dP = parallel_weierstrassPPrime(z, e1, e2, e3, nWorkers, minChunk, origSize);
    return;
end

% --- Serial core ---
dP(:) = weierPP_core(z(:).', e1(:).', e2(:).', e3(:).');


% -----------------------------------------------------------------------
function dP = weierPP_core(z, e1, e2, e3)
%WEIEPPCORE  Vectorised serial evaluation (row-vector inputs).
m    = (e2 - e3) ./ (e1 - e3);
w    = z .* sqrt(e1 - e3);
[sn, cn, dn] = ellipj(w, m);
scale = -2 .* (e1 - e3).^(3/2);
dP   = scale .* cn .* dn ./ sn.^3;
% Poles: sn -> 0 at z = 0 and at lattice points
dP(abs(sn) < eps^(1/3)) = Inf;


% -----------------------------------------------------------------------
function dP = gpu_weierstrassPPrime(z, e1, e2, e3, origSize)
%GPU_WEIERSTRASSPPRIME  GPU path via ellipj's internal GPU dispatch.
z_f  = z(:).'; e1_f = e1(:).'; e2_f = e2(:).'; e3_f = e3(:).';
m    = (e2_f - e3_f) ./ (e1_f - e3_f);
w    = z_f .* sqrt(e1_f - e3_f);
[sn, cn, dn] = ellipj(w, m);
scale = -2 .* (e1_f - e3_f).^(3/2);
dP   = scale .* cn .* dn ./ sn.^3;
dP(abs(sn) < eps^(1/3)) = Inf;
dP   = reshape(dP, origSize);


% -----------------------------------------------------------------------
function dP = parallel_weierstrassPPrime(z, e1, e2, e3, nWorkers, minChunk, origSize)
%PARALLEL_WEIERSTRASSPPRIME  Split work across parfor workers.
N       = numel(z);
z_f     = z(:).'; e1_f = e1(:).'; e2_f = e2(:).'; e3_f = e3(:).';
nChunks = min(nWorkers, ceil(N / minChunk));
chunkSz = ceil(N / nChunks);
dP_cells = cell(1, nChunks);
if exist('OCTAVE_VERSION', 'builtin')
    z_ch = cell(1,nChunks); e1_ch = cell(1,nChunks);
    e2_ch = cell(1,nChunks); e3_ch = cell(1,nChunks);
    for k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        z_ch{k}  = z_f(i1:i2);  e1_ch{k} = e1_f(i1:i2);
        e2_ch{k} = e2_f(i1:i2); e3_ch{k} = e3_f(i1:i2);
    end
    dP_cells = parcellfun(nWorkers, @par_worker, ...
        repmat({'weierstrassPPrime'}, 1, nChunks), z_ch, e1_ch, e2_ch, e3_ch, ...
        'UniformOutput', false);
else
    parfor k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        dP_cells{k} = weierstrassPPrime(z_f(i1:i2), e1_f(i1:i2), e2_f(i1:i2), e3_f(i1:i2));
    end
end
dP = reshape([dP_cells{:}], origSize);


% -----------------------------------------------------------------------
function [z, e1, e2, e3] = weierPP_broadcast(z, e1, e2, e3)
%WEIERPPBROADCAST  Expand scalar arguments to match the largest input shape.
refSz = [1 1];
for x = {z, e1, e2, e3}
    if numel(x{1}) > 1, refSz = size(x{1}); break; end
end
if isscalar(z),  z  = z( ones(refSz)); end
if isscalar(e1), e1 = e1(ones(refSz)); end
if isscalar(e2), e2 = e2(ones(refSz)); end
if isscalar(e3), e3 = e3(ones(refSz)); end
if ~(isequal(size(z),size(e1)) && isequal(size(e1),size(e2)) && isequal(size(e2),size(e3)))
    error('weierstrassPPrime: z, e1, e2, e3 must be the same size or scalar.');
end
