function [Th,H] = jacobiThetaEta(u,m,tol)
%JACOBITHETAETA evaluates Jacobi's theta and eta functions.
%   [Th, H] = JACOBITHETAETA(U,M) returns the values of the Jacobi's
%   theta and eta elliptic functions TH and H evaluated for corresponding
%   elements of argument U and parameter M.  The arrays U and M must
%   be the same size (or either can be scalar).  As currently
%   implemented, M is limited to 0 <= M <= 1.
%
%   [Th, H] = JACOBITHETAETA(U,M,TOL) computes the theta and eta
%   elliptic functions to the accuracy TOL instead of the default TOL = EPS.
%
%   Some definitions of the Jacobi elliptic functions use the modulus
%   k instead of the parameter m.  They are related by m = k^2.
%
%   Example:
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);
%       [Th, H] = jacobiThetaEta(pi/180*phi, sin(pi/180*alpha).^2);
%
%   See also
%       Standard: ELLIPKE, ELLIPJ,
%       Moiseev's package: ELLIPTIC12, ELLIPTIC12I, THETA.
%
%   ELLIPJ uses the method of the arithmetic-geometric mean
%   described in [1].
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 16-17.6.

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html
% Everyone is permitted to copy and distribute verbatim copies of this
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE.
%
% Copyright (C) 2007 by Moiseev Igor. All rights reserved.
% 34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
% For support, please reply to
%     moiseev.igor[at]gmail.com
%     Moiseev Igor,
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) || ~isreal(m)
    error('Input arguments must be real.')
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

Th = zeros(size(u));
H = Th;

if any(m(:) < 0) || any(m(:) > 1),
  error('M must be in the range 0 <= M <= 1.');
end

% GPU dispatch: move computation to GPU if enabled and available
N_el = numel(u);
if has_gpu()
    [Th,H] = gpu_jacobiThetaEta(u, m, tol);
    return;
end

% Parallel dispatch: split across workers for large inputs
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    [Th,H] = parallel_jacobiThetaEta(u, m, tol, nWorkers, minChunk);
    return;
end

m = m(:).';    % make a row vector
u = u(:).';

KK = ellipke(m);

% Theta functions from their q-series (A&S 16.27, 16.38):
%     Th(u|m) = theta_4(v, q),  H(u|m) = theta_1(v, q),  v = pi*u/(2K)
% with nome q = exp(-pi*K(1-m)/K(m)).  The series converges geometrically in
% q^(n^2) and is accurate to full double precision; the previous AGM-product
% form lost ~11 digits of the overall normalisation and needed a deliberate
% perturbation of u and m at the odd half-periods to stay finite.
q = exp(-pi .* ellipke(1-m) ./ KK);
q(~(q < 1)) = 0;                      % m == 1 (and NaN) handled below
v = pi .* u ./ (2 .* KK);

qmax = max([q(:); 0]);
if qmax > 0
    nTerms = min(1000, max(1, ceil(sqrt(log(tol) / log(qmax)))));
else
    nTerms = 1;
end

Th = ones(size(v));
H  = zeros(size(v));
for nn = 1:nTerms
    Th = Th + 2*(-1)^nn .* q.^(nn^2) .* cos(2*nn .* v);
end
for nn = 0:nTerms
    H = H + 2*(-1)^nn .* q.^((nn+0.5)^2) .* sin((2*nn+1) .* v);
end

% Special cases: m = {0, 1}
m0 = find(abs(m) < 10*eps);

if ( ~isempty(m0) )
    Th(m0) = 1;
    H(m0)  = sqrt(sqrt(m(m0))).* sin(u(m0));
end

m1 = find(abs(m-1) < 10*eps);
if ( ~isempty(m1) )
    Th(m1) = NaN;
    H(m1)  = NaN;
end


function [Th,H] = parallel_jacobiThetaEta(u, m, tol, nWorkers, minChunk)
%PARALLEL_JACOBITHETAETA  Internal helper: split work across parfor workers.
    origSize = size(u);
    N = numel(u);
    u_flat = u(:).'; m_flat = m(:).';
    nChunks = min(nWorkers, ceil(N / minChunk));
    chunkSize = ceil(N / nChunks);
    Th_c = cell(1, nChunks); H_c = cell(1, nChunks);
    if exist('OCTAVE_VERSION', 'builtin')
        u_chunks = cell(1, nChunks); m_chunks = cell(1, nChunks);
        for w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            u_chunks{w} = u_flat(i1:i2);
            m_chunks{w} = m_flat(i1:i2);
        end
        tol_c = repmat({tol}, 1, nChunks);
        tol_c = repmat({tol}, 1, nChunks);
        results = parcellfun(nWorkers, @par_worker, ...
            repmat({'jacobiThetaEta'}, 1, nChunks), u_chunks, m_chunks, tol_c, ...
            'UniformOutput', false);
        for w = 1:nChunks
            Th_c{w} = results{w}{1}; H_c{w} = results{w}{2};
        end
    else
        parfor w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            [Th_c{w}, H_c{w}] = jacobiThetaEta(u_flat(i1:i2), m_flat(i1:i2), tol);
        end
    end
    Th = reshape([Th_c{:}], origSize);
    H = reshape([H_c{:}], origSize);


function [Th,H] = gpu_jacobiThetaEta(u, m, tol)
%GPU_JACOBITHETAETA  Internal helper: q-series on gpuArray (elementwise).
%   Same series as the serial core (A&S 16.27, 16.38).  The previous GPU
%   helper still carried the retired AGM-product form and its input
%   perturbation hack, so GPU results disagreed with the CPU by up to 5e-9.
    origSize = size(u);
    Th = zeros(origSize);
    H  = zeros(origSize);
    m = m(:);  u = u(:);
    if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

    KK = ellipke(m);
    q  = exp(-pi .* ellipke(1-m) ./ KK);
    q(~(q < 1)) = 0;
    v  = pi .* u ./ (2 .* KK);
    qmax = max([q(:); 0]);
    if qmax > 0
        nTerms = min(1000, max(1, ceil(sqrt(log(tol) / log(qmax)))));
    else
        nTerms = 1;
    end
    qg = gpuArray(q);  vg = gpuArray(v);
    Thg = gpuArray(ones(size(v)));
    Hg  = gpuArray(zeros(size(v)));
    for nn = 1:nTerms
        Thg = Thg + 2*(-1)^nn .* qg.^(nn^2) .* cos(2*nn .* vg);
    end
    for nn = 0:nTerms
        Hg = Hg + 2*(-1)^nn .* qg.^((nn+0.5)^2) .* sin((2*nn+1) .* vg);
    end
    Th(:) = gather(Thg);
    H(:)  = gather(Hg);

    % Special cases: m = {0, 1}
    m0 = find(abs(m) < 10*eps);
    if ~isempty(m0), Th(m0) = 1; H(m0) = sqrt(sqrt(m(m0))) .* sin(u(m0)); end
    m1 = find(abs(m-1) < 10*eps);
    if ~isempty(m1), Th(m1) = NaN; H(m1) = NaN; end
