function Z = weierstrassZeta(z, e1, e2, e3)
%WEIERSTRASSZETA  Weierstrass zeta function ζ(z; e1, e2, e3).
%   Z = WEIERSTRASSZETA(Z, E1, E2, E3) evaluates the Weierstrass zeta
%   function, defined by  ζ'(z) = -℘(z)  and  ζ(-z) = -ζ(z)  (odd).
%
%   NOTE: this is the Weierstrass ζ, NOT the Riemann ζ.
%
%   Algorithm (closed theta form, DLMF 23.6.13 with 23.6.8):
%       ω1 = K(m)/√(e1-e3),  m = (e2-e3)/(e1-e3),  q = exp(-π·K(1-m)/K(m))
%       η1 = -π²/(12ω1) · θ1'''(0,q)/θ1'(0,q)
%       ζ(z) = η1·z/ω1 + π/(2ω1) · θ1'(v,q)/θ1(v,q),   v = πz/(2ω1)
%   The theta series converges geometrically; no quadrature is involved.
%
%   Quasi-periodicity ζ(z + 2kω1) = ζ(z) + 2kη1 is carried exactly by
%   the formula, so the function is valid for any real z.
%
%   At poles (z = 0 or lattice points): Z = Inf.
%
%   All input conventions match WEIERSTRASSP.  Parallel and GPU modes are
%   enabled via ELLIPTIC_CONFIG.  GPU acceleration cascades through the
%   internal WEIERSTRASSP calls.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions", Dover, 1965, §18.3, 18.10.
%   [2] NIST DLMF §23.6.

% Empty input -> empty output of the same shape (elementwise semantics; the
% size checks below would otherwise reject [] against a scalar).
if nargin >= 4 && (isempty(z) || isempty(e1) || isempty(e2) || isempty(e3))
    sz = size(z);
    if isempty(e1), sz = size(e1); end
    if isempty(e2), sz = size(e2); end
    if isempty(e3), sz = size(e3); end
    Z = zeros(sz);
    return;
end

if nargin < 4, error('weierstrassZeta: requires four arguments (z, e1, e2, e3).'); end
if ~isreal(z) || ~isreal(e1) || ~isreal(e2) || ~isreal(e3)
    error('weierstrassZeta: all input arguments must be real.');
end

[z, e1, e2, e3] = weierZ_broadcast(z, e1, e2, e3);
origSize = size(z);

if any(e1(:) <= e2(:)) || any(e2(:) <= e3(:))
    error('weierstrassZeta: roots must satisfy e1 > e2 > e3.');
end

Z    = zeros(origSize);
N_el = numel(z);

% GPU dispatch: serial core calls weierstrassP which auto-dispatches to GPU
if has_gpu()
    Z = gpu_weierstrassZeta(z, e1, e2, e3, origSize);
    return;
end

% Parallel dispatch
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    Z = parallel_weierstrassZeta(z, e1, e2, e3, nWorkers, minChunk, origSize);
    return;
end

% --- Serial core ---
Z(:) = weierZ_core(z(:).', e1(:).', e2(:).', e3(:).');


% -----------------------------------------------------------------------
function Z = weierZ_core(z, e1, e2, e3)
%WEIERZCORE  Vectorised serial evaluation (row-vector inputs).
%
% Closed theta form (DLMF 23.6.13 with 23.6.8):
%
%   zeta(z) = eta1*z/omega1 + pi/(2*omega1) * theta1'(v,q)/theta1(v,q)
%   eta1    = -pi^2/(12*omega1) * theta1'''(0,q)/theta1'(0,q)
%   v = pi*z/(2*omega1),  q = exp(-pi*K(1-m)/K(m)),  m = (e2-e3)/(e1-e3)
%
% No quadrature: the previous Gauss-Legendre form was ~1e-9 at best, and
% the theta series converges geometrically in q^((n+1/2)^2).  The formula
% carries the quasi-periodicity zeta(z + 2k*omega1) = zeta(z) + 2k*eta1
% exactly (theta1'(v+pi)/theta1(v+pi) is pi-periodic, the linear term does
% the rest), so no period reduction is needed either.

m_param  = (e2 - e3) ./ (e1 - e3);
mp_param = (e1 - e2) ./ (e1 - e3);    % 1-m without cancellation
one = ones(size(m_param));  zed = zeros(size(m_param));
KK  = carlsonRF(zed, mp_param, one);
KKp = carlsonRF(zed, m_param,  one);
omega1 = KK ./ sqrt(e1 - e3);
q = exp(-pi .* KKp ./ KK);
v = pi .* z ./ (2 .* omega1);

[th1, th1p, th1p0, th1ppp0] = weierZ_theta1(v, q);

eta1 = -pi^2 ./ (12 .* omega1) .* th1ppp0 ./ th1p0;
Z = eta1 .* z ./ omega1 + pi ./ (2 .* omega1) .* th1p ./ th1;

% Lattice points z = 2k*omega1 (theta1 vanishes): pole of zeta
Z(th1 == 0) = Inf;


% -----------------------------------------------------------------------
function [th1, th1p, th1p0, th1ppp0] = weierZ_theta1(v, q)
% theta1(v,q), its v-derivative, and theta1'(0), theta1'''(0).
% The common factor 2 is dropped: it cancels in every ratio used here.
qmax = max([q(:); 0]);
if qmax > 0
    nT = min(30, max(2, ceil(sqrt(abs(log(eps) ./ log(qmax))))));
else
    nT = 1;
end
th1 = zeros(size(v));  th1p = th1;
th1p0 = zeros(size(v)); th1ppp0 = th1p0;
% sin/cos of (2n+1)v by angle-addition from sin v, cos v (k*v as a double
% product rounds by eps*|k v|; see THETA_SERIES)
sk = sin(v);  ck = cos(v);  s2 = 2 .* sk .* ck;  c2 = 1 - 2 .* sk.^2;
for n = 0:nT
    qq = (-1)^n .* q.^((n+0.5)^2);
    k  = 2*n + 1;
    th1     = th1     + qq .* sk;
    th1p    = th1p    + qq .* k .* ck;
    th1p0   = th1p0   + qq .* k;
    th1ppp0 = th1ppp0 - qq .* k^3;
    [sk, ck] = deal(sk.*c2 + ck.*s2, ck.*c2 - sk.*s2);
end


% -----------------------------------------------------------------------
function Z = gpu_weierstrassZeta(z, e1, e2, e3, origSize)
%GPU_WEIERSTRASSZETA  GPU path: weierstrassP calls inside weierZ_core
%   auto-dispatch to the GPU because has_gpu() is still true.
Z = reshape(weierZ_core(z(:).', e1(:).', e2(:).', e3(:).'), origSize);


% -----------------------------------------------------------------------
function Z = parallel_weierstrassZeta(z, e1, e2, e3, nWorkers, minChunk, origSize)
%PARALLEL_WEIERSTRASSZETA  Split work across parfor workers.
N       = numel(z);
z_f     = z(:).'; e1_f = e1(:).'; e2_f = e2(:).'; e3_f = e3(:).';
nChunks = min(nWorkers, ceil(N / minChunk));
chunkSz = ceil(N / nChunks);
Z_cells = cell(1, nChunks);
if exist('OCTAVE_VERSION', 'builtin')
    z_ch = cell(1,nChunks); e1_ch = cell(1,nChunks);
    e2_ch = cell(1,nChunks); e3_ch = cell(1,nChunks);
    for k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        z_ch{k}  = z_f(i1:i2);  e1_ch{k} = e1_f(i1:i2);
        e2_ch{k} = e2_f(i1:i2); e3_ch{k} = e3_f(i1:i2);
    end
    Z_cells = parcellfun(nWorkers, @par_worker, ...
        repmat({'weierstrassZeta'}, 1, nChunks), z_ch, e1_ch, e2_ch, e3_ch, ...
        'UniformOutput', false);
else
    parfor k = 1:nChunks
        i1 = (k-1)*chunkSz + 1; i2 = min(k*chunkSz, N);
        Z_cells{k} = weierstrassZeta(z_f(i1:i2), e1_f(i1:i2), e2_f(i1:i2), e3_f(i1:i2));
    end
end
Z = reshape([Z_cells{:}], origSize);


% -----------------------------------------------------------------------
function [z, e1, e2, e3] = weierZ_broadcast(z, e1, e2, e3)
refSz = [1 1];
for x = {z, e1, e2, e3}
    if numel(x{1}) > 1, refSz = size(x{1}); break; end
end
if isscalar(z),  z  = z( ones(refSz)); end
if isscalar(e1), e1 = e1(ones(refSz)); end
if isscalar(e2), e2 = e2(ones(refSz)); end
if isscalar(e3), e3 = e3(ones(refSz)); end
if ~(isequal(size(z),size(e1)) && isequal(size(e1),size(e2)) && isequal(size(e2),size(e3)))
    error('weierstrassZeta: z, e1, e2, e3 must be the same size or scalar.');
end
