function [B, D, S] = ellipticBD(m)
%ELLIPTICBD  Complete associate elliptic integrals B(m), D(m), and S(m).
%
%   [B, D, S] = ELLIPTICBD(M) evaluates the complete associate elliptic
%   integrals:
%
%       B(m) = B(π/2|m) = ∫₀^{π/2} cos²θ / √(1−m·sin²θ) dθ
%       D(m) = D(π/2|m) = ∫₀^{π/2} sin²θ / √(1−m·sin²θ) dθ
%       S(m) = (D(m) − B(m)) / m
%
%   These are related to the standard Legendre integrals by:
%
%       K(m) = B(m) + D(m)
%       E(m) = B(m) + (1−m)·D(m) = B(m) + mc·D(m)
%       S(m) = ((2−m)·K(m) − 2·E(m)) / m²       (alternative form)
%
%   Algorithm — Carlson symmetric forms (DLMF §19.25):
%
%       B(m) = (1/2) · K(m) + (1/2) · E(m) / (1−m)   -- well-conditioned
%
%   Actually uses:
%       [K, E] = ellipke(m)
%       D(m) = (K − E) / m
%       B(m) = K − D
%       S(m) = (D − B) / m = (2D − K) / m
%
%   This avoids subtraction of nearly equal numbers via ellipke's own
%   internal cancellation-safe algorithm.
%
%   M may be a scalar or array. All elements must satisfy 0 <= m < 1.
%   At m = 0: B = D = π/4.
%   At m → 1: B → 1, D → K − 1 → ∞ (logarithmic singularity).
%
%   Parallel and GPU modes follow the library convention:
%       elliptic_config('parallel', true)
%       elliptic_config('gpu',      true)
%
%   References:
%   [1] NIST DLMF §19.2, §19.6  https://dlmf.nist.gov/19
%   [2] T. Fukushima, "Elliptic functions and elliptic integrals for
%       celestial mechanics and dynamical astronomy," (2015), §4.
%   [3] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions," Dover, 1965, §17.4.

if nargin < 1, error('ellipticBD: requires one argument m.'); end
if ~isreal(m)
    error('ellipticBD: m must be real.');
end
if any(m(:) < 0) || any(m(:) >= 1)
    error('ellipticBD: m must satisfy 0 <= m < 1.');
end

origSize = size(m);
N_el = numel(m);

% GPU dispatch
if has_gpu()
    [B, D, S] = gpu_ellipticBD(m, origSize);
    return;
end

% Parallel dispatch
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    [B, D, S] = parallel_ellipticBD(m, nWorkers, minChunk, origSize);
    return;
end

% --- Serial core ---
[B, D, S] = ellipticBD_core(m(:).', origSize);


% -----------------------------------------------------------------------
function [B, D, S] = ellipticBD_core(m, origSize)
%ELLIPTICBD_CORE  Vectorised serial evaluation (row-vector input).
[K, E] = ellipke(m);
mc = 1 - m;

% D = (K − E) / m, handle m = 0 via L'Hôpital: D(0) = π/4
D = zeros(size(m));
nz = (m ~= 0);
D(nz) = (K(nz) - E(nz)) ./ m(nz);
D(~nz) = pi / 4;

B = K - D;

% S = (D − B) / m = (2D − K) / m, handle m = 0 via L'Hôpital: S(0) = π/16
S = zeros(size(m));
S(nz) = (D(nz) - B(nz)) ./ m(nz);
S(~nz) = pi / 16;

B = reshape(B, origSize);
D = reshape(D, origSize);
S = reshape(S, origSize);


% -----------------------------------------------------------------------
function [B, D, S] = gpu_ellipticBD(m, origSize)
%GPU_ELLIPTICBD  GPU path.
[B, D, S] = ellipticBD_core(gpuArray(m(:).'), origSize);
B = gather(B);  D = gather(D);  S = gather(S);


% -----------------------------------------------------------------------
function [B, D, S] = parallel_ellipticBD(m, nWorkers, minChunk, origSize)
%PARALLEL_ELLIPTICBD  Split work across parfor workers.
N       = numel(m);
m_f     = m(:).';
nChunks = min(nWorkers, ceil(N / minChunk));
chunkSz = ceil(N / nChunks);
B_cells = cell(1, nChunks);
D_cells = cell(1, nChunks);
S_cells = cell(1, nChunks);
if exist('OCTAVE_VERSION', 'builtin')
    m_ch = cell(1, nChunks);
    for k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        m_ch{k} = m_f(i1:i2);
    end
    res = parcellfun(nWorkers, @(mv) par_worker('ellipticBD', mv), m_ch, ...
        'UniformOutput', false);
    for k = 1:nChunks
        B_cells{k} = res{k}{1};
        D_cells{k} = res{k}{2};
        S_cells{k} = res{k}{3};
    end
else
    parfor k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        [B_cells{k}, D_cells{k}, S_cells{k}] = ellipticBD(m_f(i1:i2));
    end
end
B = reshape([B_cells{:}], origSize);
D = reshape([D_cells{:}], origSize);
S = reshape([S_cells{:}], origSize);
