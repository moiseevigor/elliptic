function [B, D, J] = ellipticBDJ(phi, m, n)
%ELLIPTICBDJ  Incomplete associate elliptic integrals B(φ|m), D(φ|m), J(φ,n|m).
%
%   [B, D, J] = ELLIPTICBDJ(PHI, M, N) evaluates simultaneously:
%
%       B(φ|m) = ∫₀^φ cos²θ / √(1−m·sin²θ) dθ
%       D(φ|m) = ∫₀^φ sin²θ / √(1−m·sin²θ) dθ
%       J(φ,n|m) = ∫₀^φ sin²θ / ((1−n·sin²θ)·√(1−m·sin²θ)) dθ
%
%   [B, D] = ELLIPTICBDJ(PHI, M) computes only B and D (J = []).
%
%   Connection to standard integrals:
%
%       F(φ|m)    = B(φ|m) + D(φ|m)
%       E(φ|m)    = B(φ|m) + (1−m)·D(φ|m)
%       Π(n,φ|m)  = B(φ|m) + D(φ|m) + n·J(φ,n|m)
%
%   Computing via B,D,J separately avoids precision loss in F−E and Π−F
%   when those quantities are small.
%
%   Algorithm — Carlson symmetric forms (DLMF §19.25):
%       s = sin(φ),  c = cos(φ),  Δ = √(1−m·s²)
%
%       F  = s · R_F(c², Δ², 1)
%       D  = (s³/3) · R_D(c², Δ², 1)
%       B  = F − D
%       J  = (s³/3) · R_J(c², Δ², 1, 1−n·s²)   [if n requested]
%
%   PHI, M, N may be scalars or arrays of the same size (scalar inputs
%   are broadcast to match the largest array).  0 ≤ M < 1,  N ≠ 1.
%
%   At φ = 0: B = D = J = 0.
%   At φ = π/2: B = B(m), D = D(m) (complete associate integrals).
%
%   Parallel and GPU modes follow the library convention:
%       elliptic_config('parallel', true)
%       elliptic_config('gpu',      true)
%
%   References:
%   [1] NIST DLMF §19.25  https://dlmf.nist.gov/19.25
%   [2] T. Fukushima, "Elliptic functions and elliptic integrals for
%       celestial mechanics and dynamical astronomy," (2015), §3–4.
%   [3] B.C. Carlson, "Numerical Computation of Real or Complex Elliptic
%       Integrals," Numer. Algorithms 10 (1995), 13–26.

compute_J = (nargin >= 3);

if nargin < 2, error('ellipticBDJ: requires at least two arguments (phi, m).'); end
if ~isreal(phi) || ~isreal(m)
    error('ellipticBDJ: phi and m must be real.');
end
if compute_J && ~isreal(n)
    error('ellipticBDJ: n must be real.');
end

if compute_J
    [phi, m, n] = ellipticBDJ_broadcast3(phi, m, n);
else
    [phi, m] = ellipticBDJ_broadcast2(phi, m);
    n = [];
end
origSize = size(phi);
N_el = numel(phi);

% GPU dispatch
if has_gpu()
    [B, D, J] = gpu_ellipticBDJ(phi, m, n, compute_J, origSize);
    return;
end

% Parallel dispatch
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    [B, D, J] = parallel_ellipticBDJ(phi, m, n, compute_J, nWorkers, minChunk, origSize);
    return;
end

% --- Serial core ---
[B, D, J] = ellipticBDJ_core(phi(:).', m(:).', ...
    ifelse(compute_J, n(:).', []), compute_J, origSize);


% -----------------------------------------------------------------------
function [B, D, J] = ellipticBDJ_core(phi, m, n, compute_J, origSize)
%ELLIPTICBDJ_CORE  Vectorised serial evaluation (row-vector inputs).

s  = sin(phi);
c  = cos(phi);
d2 = 1 - m .* s.^2;    % Δ²
d  = sqrt(d2);          % Δ

% s³/3 factor
s3o3 = s.^3 ./ 3;

% R_F(c², Δ², 1) and R_D(c², Δ², 1)
RF = carlsonRF(c.^2, d2, ones(size(phi)));
RD = carlsonRD(c.^2, d2, ones(size(phi)));

F_val = s .* RF;          % = F(φ|m)
D_val = s3o3 .* RD;       % = D(φ|m)
B_val = F_val - D_val;    % = B(φ|m)

% Poles / zeros at s=0
B_val(s == 0) = 0;
D_val(s == 0) = 0;

B = reshape(B_val, origSize);
D = reshape(D_val, origSize);

if compute_J
    % 1 − n·s² (denominator parameter for R_J)
    p = 1 - n .* s.^2;
    RJ = carlsonRJ(c.^2, d2, ones(size(phi)), p);
    J_val = s3o3 .* RJ;
    J_val(s == 0) = 0;
    J = reshape(J_val, origSize);
else
    J = [];
end


% -----------------------------------------------------------------------
function [B, D, J] = gpu_ellipticBDJ(phi, m, n, compute_J, origSize)
[B, D, J] = ellipticBDJ_core(gpuArray(phi(:).'), gpuArray(m(:).'), ...
    ifelse(compute_J, gpuArray(n(:).'), []), compute_J, origSize);
B = gather(B); D = gather(D);
if compute_J, J = gather(J); end


% -----------------------------------------------------------------------
function [B, D, J] = parallel_ellipticBDJ(phi, m, n, compute_J, nWorkers, minChunk, origSize)
N       = numel(phi);
phi_f   = phi(:).'; m_f = m(:).';
if compute_J, n_f = n(:).'; else n_f = []; end
nChunks = min(nWorkers, ceil(N / minChunk));
chunkSz = ceil(N / nChunks);
B_cells = cell(1,nChunks); D_cells = cell(1,nChunks); J_cells = cell(1,nChunks);
if exist('OCTAVE_VERSION', 'builtin')
    phi_ch=cell(1,nChunks); m_ch=cell(1,nChunks); n_ch=cell(1,nChunks);
    for k=1:nChunks
        i1=(k-1)*chunkSz+1; i2=min(k*chunkSz,N);
        phi_ch{k}=phi_f(i1:i2); m_ch{k}=m_f(i1:i2);
        if compute_J, n_ch{k}=n_f(i1:i2); end
    end
    if compute_J
        res=parcellfun(nWorkers,@(p,mv,nv) par_worker('ellipticBDJ',p,mv,nv),...
            phi_ch,m_ch,n_ch,'UniformOutput',false);
    else
        res=parcellfun(nWorkers,@(p,mv) par_worker('ellipticBDJ',p,mv),...
            phi_ch,m_ch,'UniformOutput',false);
    end
    for k=1:nChunks
        B_cells{k}=res{k}{1}; D_cells{k}=res{k}{2};
        if compute_J, J_cells{k}=res{k}{3}; end
    end
else
    parfor k=1:nChunks
        i1=(k-1)*chunkSz+1; i2=min(k*chunkSz,N);
        if compute_J
            [B_cells{k},D_cells{k},J_cells{k}]=ellipticBDJ(phi_f(i1:i2),m_f(i1:i2),n_f(i1:i2));
        else
            [B_cells{k},D_cells{k}]=ellipticBDJ(phi_f(i1:i2),m_f(i1:i2));
        end
    end
end
B = reshape([B_cells{:}], origSize);
D = reshape([D_cells{:}], origSize);
if compute_J, J = reshape([J_cells{:}], origSize); else J = []; end


% -----------------------------------------------------------------------
function [phi, m] = ellipticBDJ_broadcast2(phi, m)
refSz = [1 1];
for v = {phi, m}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(phi), phi = phi(ones(refSz)); end
if isscalar(m),   m   = m(ones(refSz));   end
if ~isequal(size(phi), size(m))
    error('ellipticBDJ: phi and m must be the same size or scalar.');
end


% -----------------------------------------------------------------------
function [phi, m, n] = ellipticBDJ_broadcast3(phi, m, n)
refSz = [1 1];
for v = {phi, m, n}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(phi), phi = phi(ones(refSz)); end
if isscalar(m),   m   = m(ones(refSz));   end
if isscalar(n),   n   = n(ones(refSz));   end
if ~(isequal(size(phi),size(m)) && isequal(size(m),size(n)))
    error('ellipticBDJ: phi, m, n must be the same size or scalar.');
end


% -----------------------------------------------------------------------
function y = ifelse(cond, a, b)
if cond, y = a; else y = b; end
