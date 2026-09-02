function Pi = elliptic3(u,m,c);
% ELLIPTIC3 evaluates incomplete elliptic integral of the third kind.
%   Pi = ELLIPTIC3(U,M,C) where U is a phase in radians, 0<M<1 is
%   the module and 0<C<1 is a parameter.
%
%   ELLIPTIC3 uses 20-node Gauss-Legendre quadrature away from endpoint
%   singularities and Carlson RF/RJ forms near a pole.  The hybrid avoids
%   the quadrature's loss of precision as C or M approaches one while
%   retaining its vectorised fast path for regular inputs.
%
%   Pi(u,m,c) = int(1/((1 - c*sin(t)^2)*sqrt(1 - m*sin(t)^2)), t=0..u)
%
%   Tables generating code ([1], pp. 625-626):
%	    [phi,alpha,c] = meshgrid(0:15:90, 0:15:90, 0:0.1:1);
%   	Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 17.7.
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989.
%   [3] S. Zhang, J. Jin "Computation of Special Functions" (Wiley, 1996).

%   For support, please reply to
%       moiseev.igor[at]gmail.com
%       Moiseev Igor,
%       34106, SISSA, via Beirut n. 2-4,  Trieste, Italy

if nargin<3, error('Not enough input arguments.'); end
if ~isreal(u) || ~isreal(m) || ~isreal(c)
    error('Input arguments must be real.')
end
if any(m < 0) || any(m > 1) || any(c > 1),
  error('M must be in the range [0, 1] and C <= 1.');
end
% Reduce the phase to [0, pi/2] using oddness and the quasi-period
% (the integrand is pi-periodic and even about every multiple of pi/2):
%   Pi(-u|m,c)     = -Pi(u|m,c)
%   Pi(u+k*pi|m,c) =  Pi(u|m,c) + 2k*Pi(pi/2|m,c)
%   Pi(pi-u|m,c)   =  2*Pi(pi/2|m,c) - Pi(u|m,c)
% The Gauss-Legendre rule below is only accurate on [0, pi/2]; phases
% outside that range used to be rejected with an error.
if any(u(:) < 0) || any(u(:) > pi/2)
    signU = sign(u);  ua = abs(u);
    k_per = floor(ua ./ pi);
    r     = sub_kpi(ua, k_per);                                     % in [0, pi), error eps*|r| (see SUB_KPI)
    refl  = r > pi/2;
    ur    = r;  ur(refl) = pi - r(refl); % in [0, pi/2]
    Pred  = elliptic3(ur, m, c);
    % Complete-integral correction only where a half-period or reflection
    % actually applies: when the complete integral is a pole (m = 1 or
    % c = 1), an unconditional 2*k_per*Pcpl forms 0*Inf = NaN for phases
    % that never cross it (e.g. plain negative amplitudes).
    corr = zeros(size(ur));
    idx  = (k_per > 0) | refl;
    if any(idx(:))
        % Exact complete integral Pi(c|m) = R_F(0,1-m,1) + (c/3) R_J(0,1-m,1,1-c)
        % (DLMF 19.25.?) -- NOT elliptic3(pi/2,...): cos(double(pi/2)) = 6e-17 is
        % not 0, and near m = 1 the sliver between double(pi/2) and pi/2 is
        % 6e-17/(sqrt(1-m)(1-c)) ~ 2e-7 (relative 3e-10 in the reflected value).
        Pcpl = carlsonRF(zeros(size(u)), 1 - m, ones(size(u))) + ...
               c ./ 3 .* carlsonRJ(zeros(size(u)), 1 - m, ones(size(u)), 1 - c);
        corr(idx) = 2 .* k_per(idx) .* Pcpl(idx) + 2 .* refl(idx) .* Pcpl(idx);
    end
    Pi = signU .* (corr + (1 - 2 .* refl) .* Pred);
    return;
end

[mm,nm] = size(m);
[mu,nu] = size(u);
% Broadcast scalars to the largest input (the old order expanded c from the
% still-scalar u before u itself was expanded from m, so a scalar phase with a
% parameter vector was rejected as 'must be the same size').
sz = size(u);
if numel(m) > 1, sz = size(m); elseif numel(c) > 1, sz = size(c); end
if isempty(u) || isempty(m) || isempty(c), Pi = zeros(0, 0); if isempty(u), Pi = zeros(size(u)); elseif isempty(m), Pi = zeros(size(m)); else, Pi = zeros(size(c)); end; return; end
if length(m)==1, m = m(ones(sz)); end
if length(c)==1, c = c(ones(sz)); end
if length(u)==1, u = u(ones(sz)); end
if ~isequal(size(m), size(c), size(u)),
        error('U, M and C must be the same size.');
end

Pi = zeros(size(u));

% GPU dispatch: the legacy fixed quadrature is retained only away from the
% endpoint singularities.  Near a pole it loses several digits, so fall back
% to the Carlson CPU core below.
N_el = numel(u);
gpu_regular = all((1 - c(:).*sin(u(:)).^2) >= 0.25) && ...
              all((1 - m(:).*sin(u(:)).^2) >= 0.25);
if has_gpu() && gpu_regular
    Pi = gpu_elliptic3(u, m, c);
    return;
end

% Parallel dispatch: split across workers for large inputs
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    Pi = parallel_elliptic3(u, m, c, nWorkers, minChunk);
    return;
end

m = m(:).';    % make a row vector
u = u(:).';
c = c(:).';

I = find(u==pi/2 & m==1 | u==pi/2 & c==1);

% Hybrid evaluator.  The 20-node rule is full precision while both endpoint
% denominators stay >= 0.25; nearer a pole, switch only those elements to the
% Carlson form (DLMF 19.25.14).  This retains the vectorised fast path without
% the previous seven-digit loss as c approached 1.
s = sin(u);
s2 = s.^2;
co = cos(u);
% (1-m) + m cos^2 and (1-c) + c cos^2: no cancellation near the endpoint
% poles (Pi(pi/2-1e-6 | m, c=1) was off by 3e-5).
d2 = (1 - m) + m.*co.^2;
p = (1 - c) + c.*co.^2;
% c < 0 (allowed, as in the Python port): the integrand 1/(1 + |c| sin^2)
% narrows as |c| grows and the 20-node rule loses digits, so use Carlson.
danger = (d2 < 0.25) | (p < 0.25) | (c < 0);
P = zeros(size(u));

regular = find(~danger);
if ~isempty(regular)
    t = [0.9931285991850949, 0.9639719272779138, ...
         0.9122344282513259, 0.8391169718222188, ...
         0.7463319064601508, 0.6360536807265150, ...
         0.5108670019508271, 0.3737060887154195, ...
         0.2277858511416451, 0.07652652113349734];
    w = [0.01761400713915212, 0.04060142980038694, ...
         0.06267204833410907, 0.08327674157670475, ...
         0.1019301198172404,  0.1181945319615184, ...
         0.1316886384491766,  0.1420961093183820, ...
         0.1491729864726037,  0.1527533871307258];
    ur = u(regular);
    mr = m(regular);
    cr = c(regular);
    Pr = zeros(size(ur));
    for jj = 1:10
        c0 = ur .* t(jj) ./ 2;
        Pr = Pr + w(jj) .* (g(ur./2+c0, mr, cr) + g(ur./2-c0, mr, cr));
    end
    P(regular) = ur ./ 2 .* Pr;
end

% Keep the eager Carlson evaluation finite at endpoint poles; those outputs
% are replaced by Inf below.
near = find(danger);
if ~isempty(near)
    P(near) = elliptic3_carlson(u(near), m(near), c(near));
end
P(s == 0) = 0;
Pi(:) = P;

% special values u==pi/2 & m==1 | u==pi/2 & c==1
Pi(I) = inf;
return;


function P = elliptic3_carlson(u, m, c)
%ELLIPTIC3_CARLSON  Pi(u|m,c) by DLMF 19.25.14 for 0 <= u <= pi/2 (row inputs).
%   Used for the elements the 20-node rule cannot resolve (denominators
%   below 0.25 at the endpoint, or c < 0) by BOTH the serial core and the
%   GPU path: the OpenCL kernel is the quadrature only, and on the L4 it
%   returned Pi(1|0.5,-100) 3.8e-9 off because it had no such fallback.
s  = sin(u);  co = cos(u);
d2 = (1 - m) + m.*co.^2;
p  = (1 - c) + c.*co.^2;
endpoint = (u == pi/2 & (m == 1 | c == 1));   % keep the eager evaluation finite; caller sets Inf
c(endpoint) = 0;  d2(endpoint) = 1;  p(endpoint) = 1;
RF = carlsonRF(co.^2, d2, ones(size(u)));
RJ = carlsonRJ(co.^2, d2, ones(size(u)), p);
P  = s.*RF + c.*s.^3.*RJ./3;
P(s == 0) = 0;


function g = g(u,m,c)
%  g = 1/((1 - c*sin(u)^2)*sqrt(1 - m*sin(u)^2));

 cs2 = cos(u).^2;
 g = 1./(((1 - c) + c.*cs2).*sqrt((1 - m) + m.*cs2));
return;


function Pi = parallel_elliptic3(u, m, c, nWorkers, minChunk)
%PARALLEL_ELLIPTIC3  Internal helper: split work across parfor workers.
    origSize = size(u);
    N = numel(u);
    u_flat = u(:).'; m_flat = m(:).'; c_flat = c(:).';
    nChunks = min(nWorkers, ceil(N / minChunk));
    chunkSize = ceil(N / nChunks);
    Pi_cells = cell(1, nChunks);
    if exist('OCTAVE_VERSION', 'builtin')
        u_chunks = cell(1, nChunks);
        m_chunks = cell(1, nChunks);
        c_chunks = cell(1, nChunks);
        for w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            u_chunks{w} = u_flat(i1:i2);
            m_chunks{w} = m_flat(i1:i2);
            c_chunks{w} = c_flat(i1:i2);
        end
        Pi_cells = parcellfun(nWorkers, @par_worker, ...
            repmat({'elliptic3'}, 1, nChunks), u_chunks, m_chunks, c_chunks, ...
            'UniformOutput', false);
    else
        parfor w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            Pi_cells{w} = elliptic3(u_flat(i1:i2), m_flat(i1:i2), c_flat(i1:i2));
        end
    end
    Pi = reshape([Pi_cells{:}], origSize);


function Pi = gpu_elliptic3(u, m, c)
%GPU_ELLIPTIC3  Internal helper: compute elliptic3 using gpuArray (MATLAB only).
%   The Gauss-Legendre integrand is elementwise — trivially GPU-compatible.
    origSize = size(u);
    I_inf = find(u(:).' == pi/2 & m(:).' == 1 | u(:).' == pi/2 & c(:).' == 1);

    uu = u(:).';  mm = m(:).';  cc = c(:).';
    % Same hybrid as the serial core: the kernel is the 20-node rule, which is
    % full precision only while both endpoint denominators stay >= 0.25 and
    % c >= 0; the rest goes through the Carlson form on the host.
    co2 = cos(uu).^2;
    danger = ((1 - mm) + mm.*co2 < 0.25) | ((1 - cc) + cc.*co2 < 0.25) | (cc < 0) | isnan(uu) | isnan(mm) | isnan(cc);
    Pi = zeros(origSize);
    if any(danger)
        Pi(danger) = elliptic3_carlson(uu(danger), mm(danger), cc(danger));
        Pi(isnan(uu) | isnan(mm) | isnan(cc)) = NaN;
    end
    if ~any(~danger)
        Pi(I_inf) = inf;
        return;
    end
    reg = ~danger;
    u_g = gpuArray(uu(reg));
    m_g = gpuArray(mm(reg));
    c_g = gpuArray(cc(reg));

    t = [ 0.9931285991850949,  0.9639719272779138, ...
          0.9122344282513259,  0.8391169718222188, ...
          0.7463319064601508,  0.6360536807265150, ...
          0.5108670019508271,  0.3737060887154195, ...
          0.2277858511416451,  0.07652652113349734 ];
    w = [ 0.01761400713915212, 0.04060142980038694, ...
          0.06267204833410907, 0.08327674157670475, ...
          0.1019301198172404,  0.1181945319615184,  ...
          0.1316886384491766,  0.1420961093183820,  ...
          0.1491729864726037,  0.1527533871307258   ];

    P = gpuArray(zeros(size(u_g)));
    for ii = 1:10
        c0 = u_g .* t(ii) / 2;
        P  = P + w(ii) .* (g_gpu(u_g/2 + c0, m_g, c_g) + g_gpu(u_g/2 - c0, m_g, c_g));
    end
    P = u_g/2 .* P;

    Pi(reg) = gather(P);
    Pi(I_inf) = inf;


function gv = g_gpu(u, m, c)
    cs2 = cos(u).^2;
    gv  = 1 ./ (((1 - c) + c.*cs2) .* sqrt((1 - m) + m.*cs2));
