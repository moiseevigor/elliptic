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
period_condition = u./KK/2-floor(u./KK/2);

I_odd = find( abs(m-1) > 10*eps & abs(m) > 10*eps & abs(period_condition - 0.5) < 10*eps );
% here we cheat, add some disturbance to avoid the AGM algorithm divergence
% when the inputs are the ratios of complete elliptic integrals
if ( ~isempty(I_odd) )
    u(I_odd) = u(I_odd) + 100000*eps;
    m(I_odd) = m(I_odd) + 10000*eps;
end

I = find( abs(m-1) > 10*eps & ...
          abs(m) > 10*eps ...
        );
%         abs(period_condition - 0.5) > 10*eps ...                % odd period
%         abs(period_condition) > 10*eps ...                      % even period

if ~isempty(I)
    % Use standard uniquetol for numerical precision issues
    % This is the recommended MATLAB approach since R2015a
    m_vals = m(I);
    tol_unique = 1e-11;

    [mu, ~, K] = uniquetol_compat(m_vals, tol_unique);
    K = K(:).';

    % pre-allocate space and augment if needed
	chunk = 7;
	a = zeros(chunk,length(mu));
	c = a;
	b = a;
	a(1,:) = ones(1,length(mu));
	c(1,:) = sqrt(mu);
	b(1,:) = sqrt(1-mu);
	n = zeros(1,length(mu));
	i = 1;

    % Arithmetic-Geometric Mean of A, B and C
    while any(abs(c(i,:)) > tol)
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,length(mu))];
          b = [b; zeros(2,length(mu))];
          c = [c; zeros(2,length(mu))];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        mask = (abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol);
        n(mask) = i-1;
	end

    mmax = length(I);
	phin = zeros(1,mmax);
    prodth = ones(i,mmax);

    % Calculate phin
    phin(:) = (2 .^ n(K)).*a(i,K).*u(I);
    phin_pred = phin;
	while i > 1
        i = i - 1;
        mask = n(K) >= i;
        if any(mask)
          phin(mask) = 0.5*(asin(c(i+1,K(mask)).*sin(phin(mask))./a(i+1,K(mask))) + phin(mask));
          prodth(i,mask) = ( sec(2*phin(mask)-phin_pred(mask)) ).^(1/2^(i+1));
          if (i > 1)
              phin_pred = phin;
          end
        end
    end

    th_save = sqrt(2*sqrt(1-m(I)).* KK(I)/pi.* cos(phin_pred - phin)./cos(phin) ).* prod(prodth,1);
    Th(I) = th_save;
    H(I) = sqrt(sqrt(m(I))).* sin(phin).* th_save;
end

% special values of u = (2n+1)*KK, odd periods
% I_odd = find( abs(m-1) > 10*eps & abs(m) > 10*eps & abs(period_condition - 0.5) < 10*eps );
% if ( ~isempty(I_odd) )
%     Th(I_odd) = 1+2*(1./(1-exp(-pi*ellipke(1-m(I_odd))./KK(I_odd)))-1);
%     Th(I_odd) = sqrt(KK(I_odd)./ellipke(1-m(I_odd)));
%     H(I_odd)  = (-1).^(floor(u(I_odd)./KK(I_odd)/2)).* sqrt(sqrt(m(I_odd))).* Th(I_odd);
% end

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
%GPU_JACOBITHETAETA  Internal helper: compute jacobiThetaEta using gpuArray.
%   Compatible with both MATLAB gpuArray and Octave ocl package.
    origSize = size(u);
    Th = zeros(origSize);
    H  = zeros(origSize);

    m = m(:);
    u = u(:);

    if any(m(:) < 0) || any(m(:) > 1), error('M must be in the range 0 <= M <= 1.'); end

    KK = ellipke(m);
    period_condition = u./KK/2 - floor(u./KK/2);

    I_odd = find(abs(m-1) > 10*eps & abs(m) > 10*eps & abs(period_condition - 0.5) < 10*eps);
    if ~isempty(I_odd)
        u(I_odd) = u(I_odd) + 100000*eps;
        m(I_odd) = m(I_odd) + 10000*eps;
    end

    I = find(abs(m-1) > 10*eps & abs(m) > 10*eps);
    if ~isempty(I)
        mmax = length(I);
        mu   = m(I);

        % Transposed layout: rows=elements, cols=iterations (OCL-friendly)
        MAX_ITER = 12;
        a = gpuArray(zeros(mmax, MAX_ITER));
        b = gpuArray(zeros(mmax, MAX_ITER));
        c = gpuArray(zeros(mmax, MAX_ITER));
        a(:,1) = gpuArray(ones(mmax,1));
        c(:,1) = gpuArray(sqrt(mu));
        b(:,1) = gpuArray(sqrt(1 - mu));
        n = zeros(mmax, 1);
        ii = 1;
        while any(gather(abs(c(:,ii))) > tol)
            ii = ii + 1;
            a(:,ii) = 0.5 * (a(:,ii-1) + b(:,ii-1));
            b(:,ii) = sqrt(a(:,ii-1) .* b(:,ii-1));
            c(:,ii) = 0.5 * (a(:,ii-1) - b(:,ii-1));
            mask = logical(gather((abs(c(:,ii)) <= tol) & (abs(c(:,ii-1)) > tol)));
            n(mask) = ii - 1;
        end

        % Ascending Landen back-substitution with multiplicative masking
        phin      = gpuArray((2 .^ n) .* gather(a(:,ii)) .* u(I));
        phin_pred = phin;
        prodth    = gpuArray(ones(mmax, MAX_ITER));
        for jj = ii-1:-1:1
            active = gpuArray(double(n >= jj));
            phin_new = 0.5*(asin(c(:,jj+1).*sin(phin)./a(:,jj+1)) + phin);
            phin_upd = phin + active .* (phin_new - phin);
            prodth(:,jj) = 1 + active .* ((sec(2*phin_upd - phin_pred)).^(1/2^(jj+1)) - 1);
            if jj > 1, phin_pred = phin_upd; end
            phin = phin_upd;
        end

        th_save = sqrt(2*sqrt(1 - m(I)) .* KK(I)/pi .* ...
            gather(cos(phin_pred - phin) ./ cos(phin))) .* gather(prod(prodth, 2));
        Th(I) = th_save;
        H(I)  = sqrt(sqrt(m(I))) .* gather(sin(phin)) .* th_save;
    end

    % Special cases: m = {0, 1}
    m0 = find(abs(m) < 10*eps);
    if ~isempty(m0), Th(m0) = 1; H(m0) = sqrt(sqrt(m(m0))) .* sin(u(m0)); end

    m1 = find(abs(m-1) < 10*eps);
    if ~isempty(m1), Th(m1) = NaN; H(m1) = NaN; end