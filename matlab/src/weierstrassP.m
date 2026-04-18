function P = weierstrassP(z, e1, e2, e3)
%WEIERSTRASSP  Weierstrass elliptic function ℘(z; e1, e2, e3).
%   P = WEIERSTRASSP(Z, E1, E2, E3) evaluates the Weierstrass P-function
%   at argument Z for the lattice whose half-period values are E1 > E2 > E3
%   (the three roots of  4t^3 - g2*t - g3 = 0).
%
%   Formula (A&S 18.9.1):
%       m  = (e2 - e3) / (e1 - e3)
%       w  = z * sqrt(e1 - e3)
%       P  = e3 + (e1 - e3) / sn^2(w, m)
%
%   Z, E1, E2, E3 must be real.  Any combination of scalar / same-size
%   arrays is accepted; scalars are broadcast to match the largest input.
%   E1+E2+E3 = 0 is recommended (standard lattice normalisation, A&S 18.1.2)
%   but not enforced.
%
%   At lattice poles (z = 0 or any period), P = Inf.
%
%   Parallel and GPU acceleration follow the library convention:
%       elliptic_config('parallel', true)   % multi-core via parcellfun/parfor
%       elliptic_config('gpu',      true)   % GPU via gpuArray (needs ellipj GPU)
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions", Dover, 1965, §18.9.
%   [2] NIST DLMF §23.6.

if nargin < 4, error('weierstrassP: requires four arguments (z, e1, e2, e3).'); end
if ~isreal(z) || ~isreal(e1) || ~isreal(e2) || ~isreal(e3)
    error('weierstrassP: all input arguments must be real.');
end

% Scalar broadcasting: expand any singleton to match the largest input
[z, e1, e2, e3] = weierP_broadcast(z, e1, e2, e3);
origSize = size(z);

if any(e1(:) <= e2(:)) || any(e2(:) <= e3(:))
    error('weierstrassP: roots must satisfy e1 > e2 > e3.');
end

P    = zeros(origSize);
N_el = numel(z);

% GPU dispatch: move computation to GPU if enabled and available
if has_gpu()
    P = gpu_weierstrassP(z, e1, e2, e3, origSize);
    return;
end

% Parallel dispatch: split across workers for large inputs
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    P = parallel_weierstrassP(z, e1, e2, e3, nWorkers, minChunk, origSize);
    return;
end

% --- Serial core ---
P(:) = weierP_core(z(:).', e1(:).', e2(:).', e3(:).');


% -----------------------------------------------------------------------
function P = weierP_core(z, e1, e2, e3)
%WEIEP_CORE  Vectorised serial evaluation (row-vector inputs).
m   = (e2 - e3) ./ (e1 - e3);
w   = z .* sqrt(e1 - e3);
[sn, ~, ~] = ellipj(w, m);
P   = e3 + (e1 - e3) ./ sn.^2;
% Poles: sn -> 0 at z = 0 and at lattice points
P(abs(sn) < eps^(1/3)) = Inf;


% -----------------------------------------------------------------------
function P = gpu_weierstrassP(z, e1, e2, e3, origSize)
%GPU_WEIERSTRASSP  GPU path: ellipj handles its own GPU dispatch internally.
z_f  = z(:).'; e1_f = e1(:).'; e2_f = e2(:).'; e3_f = e3(:).';
m    = (e2_f - e3_f) ./ (e1_f - e3_f);
w    = z_f .* sqrt(e1_f - e3_f);
% ellipj sees has_gpu()=true and dispatches to gpu_ellipj automatically
[sn, ~, ~] = ellipj(w, m);
P    = e3_f + (e1_f - e3_f) ./ sn.^2;
P(abs(sn) < eps^(1/3)) = Inf;
P    = reshape(P, origSize);


% -----------------------------------------------------------------------
function P = parallel_weierstrassP(z, e1, e2, e3, nWorkers, minChunk, origSize)
%PARALLEL_WEIERSTRASSP  Split work across parfor workers.
N       = numel(z);
z_f     = z(:).'; e1_f = e1(:).'; e2_f = e2(:).'; e3_f = e3(:).';
nChunks = min(nWorkers, ceil(N / minChunk));
chunkSz = ceil(N / nChunks);
P_cells = cell(1, nChunks);
if exist('OCTAVE_VERSION', 'builtin')
    z_ch = cell(1,nChunks); e1_ch = cell(1,nChunks);
    e2_ch = cell(1,nChunks); e3_ch = cell(1,nChunks);
    for k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        z_ch{k}  = z_f(i1:i2);  e1_ch{k} = e1_f(i1:i2);
        e2_ch{k} = e2_f(i1:i2); e3_ch{k} = e3_f(i1:i2);
    end
    P_cells = parcellfun(nWorkers, @par_worker, ...
        repmat({'weierstrassP'}, 1, nChunks), z_ch, e1_ch, e2_ch, e3_ch, ...
        'UniformOutput', false);
else
    parfor k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        P_cells{k} = weierstrassP(z_f(i1:i2), e1_f(i1:i2), e2_f(i1:i2), e3_f(i1:i2));
    end
end
P = reshape([P_cells{:}], origSize);


% -----------------------------------------------------------------------
function [z, e1, e2, e3] = weierP_broadcast(z, e1, e2, e3)
%WEIEP_BROADCAST  Expand scalar arguments to match the largest input shape.
inputs = {z, e1, e2, e3};
% Find the reference (non-scalar) size
refSz = [1 1];
for k = 1:4
    if numel(inputs{k}) > 1, refSz = size(inputs{k}); break; end
end
if isscalar(z),  z  = z( ones(refSz));  end
if isscalar(e1), e1 = e1(ones(refSz)); end
if isscalar(e2), e2 = e2(ones(refSz)); end
if isscalar(e3), e3 = e3(ones(refSz)); end
if ~(isequal(size(z),size(e1)) && isequal(size(e1),size(e2)) && isequal(size(e2),size(e3)))
    error('weierstrassP: z, e1, e2, e3 must be the same size or scalar.');
end
