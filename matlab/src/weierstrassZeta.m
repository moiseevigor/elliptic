function Z = weierstrassZeta(z, e1, e2, e3)
%WEIERSTRASSZETA  Weierstrass zeta function ζ(z; e1, e2, e3).
%   Z = WEIERSTRASSZETA(Z, E1, E2, E3) evaluates the Weierstrass zeta
%   function, defined by  ζ'(z) = -℘(z)  and  ζ(-z) = -ζ(z)  (odd).
%
%   NOTE: this is the Weierstrass ζ, NOT the Riemann ζ.
%
%   Algorithm:
%     1. Compute the real half-period  ω1 = K(m)/√(e1-e3),  m=(e2-e3)/(e1-e3).
%     2. Compute the quasi-period  η1 = ζ(ω1)  via a regularised integral
%        (A&S 18.10.1):
%            η1 = 1/ω1 + ∫_0^{ω1} [-℘(t) + 1/t²] dt
%        The integrand equals O(t²) near t=0 (Laurent cancellation), so it
%        is evaluated with a 10-point Gauss-Legendre rule on [ε₀, ω1] plus
%        an analytic correction  -g2*ε₀³/60 - g3*ε₀⁵/140  on [0, ε₀]
%        (ε₀ = 0.1·ω1 avoids catastrophic cancellation).
%     3. Integrate:  ζ(z) = η1 + ∫_{ω1}^z [-℘(t)] dt  (GL quadrature).
%
%   Quasi-periodicity: ζ(z + 2ω1) = ζ(z) + 2η1.  The function is valid
%   for |z| < ~4ω1; accuracy degrades near lattice poles.
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
N = numel(z);

% Half-period and modulus
m_param = (e2 - e3) ./ (e1 - e3);
[KK, ~] = ellipke(m_param);
omega1  = KK ./ sqrt(e1 - e3);

% Lattice invariants needed for the small-t analytic correction
g2 = -4 .* (e1.*e2 + e1.*e3 + e2.*e3);
g3 =  4 .* e1 .* e2 .* e3;

% η1 = ζ(ω1) via regularised integral on [ε₀, ω1] + analytic patch on [0, ε₀]
eps0   = 0.1 .* omega1;                     % cutoff: avoids cancellation
% Analytic correction for ∫_0^{ε₀} [-℘(t)+1/t²] dt:
%   -℘(t)+1/t² = -g2*t²/20 - g3*t⁴/28 + ...  (A&S 18.5.5)
%   integral   = -g2*ε₀³/60 - g3*ε₀⁵/140
anal   = -g2 .* eps0.^3 / 60 - g3 .* eps0.^5 / 140;
% GL integral on [ε₀, ω1] of  f(t) = -℘(t) + 1/t²
eta1   = 1 ./ omega1 + anal + weierZ_gl_integral(@feta, eps0, omega1, e1, e2, e3);

% Period reduction: bring z to (-ω1, ω1] so the integration path ω1→z_red
% stays within (-ω1, ω1] and never crosses a pole of ℘ at 0 or ±2ω1, etc.
%   ζ(z) = ζ(z_red) + k·2η1   where  z_red = z - k·2ω1
% Use floor with a half-period offset to map (-ω1, ω1] (avoids round(0.5)=1
% pushing ω1 itself to the wrong strip).
k      = floor((z + omega1 .* (1 - eps)) ./ (2 .* omega1));
z_red  = z - k .* (2 .* omega1);

% ζ(z_red) = η1 + ∫_{ω1}^{z_red} [-℘(t)] dt, but only when z_red > 0.
% For z_red < 0 the path from ω1 to z_red crosses the pole of ℘ at t=0.
% Use odd symmetry ζ(-z) = -ζ(z) instead: evaluate at |z_red| and negate.
neg    = z_red < 0;
az_red = abs(z_red);
int_val = weierZ_gl_integral(@fzeta, omega1, az_red, e1, e2, e3);
Z_red = eta1 + int_val;          % ζ(|z_red|) for all elements
Z_red(neg) = -Z_red(neg);        % ζ(z_red) = -ζ(-z_red) for negative z_red
Z = Z_red + k .* (2 .* eta1);

% Poles: z_red ≈ 0 means z is at a lattice point
Z(az_red < eps^(1/3)) = Inf;


% -----------------------------------------------------------------------
function v = feta(t, e1, e2, e3)
% Regularised integrand for η1: -℘(t) + 1/t²  (regular at t=0, O(t²))
P = weierstrassP(t, e1, e2, e3);
v = -P + 1 ./ t.^2;
% At extremely small t the lattice pole of P makes ℘ ≈ 1/t², so the
% difference can underflow; replace with analytic O(t²) value if needed.
g2 = -4 .* (e1.*e2 + e1.*e3 + e2.*e3);
g3 =  4 .* e1 .* e2 .* e3;
small = abs(t) < 1e-4;
if any(small)
    v(small) = (-g2(small) .* t(small).^2 / 20 - g3(small) .* t(small).^4 / 28);
end


% -----------------------------------------------------------------------
function v = fzeta(t, e1, e2, e3)
% Integrand for ζ(z): -℘(t)
v = -weierstrassP(t, e1, e2, e3);
% Replace Inf (at poles) with 0 so GL sum doesn't blow up;
% in practice t never hits an exact lattice point in the GL grid.
v(~isfinite(v)) = 0;


% -----------------------------------------------------------------------
function I = weierZ_gl_integral(f_handle, a, b, e1, e2, e3)
%WEIERZ_GL_INTEGRAL  10-pt GL quadrature of f(t,e1,e2,e3) from a to b (arrays).
%   All of a, b, e1, e2, e3 are row vectors of the same length N.
%   Returns a row vector I of length N.
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

mid = (a + b) ./ 2;
hw  = (b - a) ./ 2;
I   = zeros(size(a));
for k = 1:10
    tp = mid + hw .* t_gl(k);
    tm = mid - hw .* t_gl(k);
    I  = I + w_gl(k) .* (f_handle(tp, e1, e2, e3) + f_handle(tm, e1, e2, e3));
end
I = I .* hw;


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
