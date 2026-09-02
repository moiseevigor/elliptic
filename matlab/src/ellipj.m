function [sn,cn,dn,am] = ellipj(u,m,tol)
%ELLIPJ Jacobi elliptic functions and Jacobi's amplitude.
%   [Sn,Cn,Dn,Am] = ELLIPJ(U,M) returns the values of the Jacobi
%   elliptic functions SN, CN, DN and AM evaluated for corresponding
%   elements of argument U and parameter M.  The arrays U and M must
%   be the same size (or either can be scalar).  As currently
%   implemented, M is limited to 0 <= M <= 1.
%
%   [Sn,Cn,Dn,Am] = ELLIPJ(U,M,TOL) computes the elliptic functions to
%   the accuracy TOL instead of the default TOL = EPS.
%
%   Accuracy limit for large arguments: the phase is reduced modulo 2K in
%   double precision, so the residual phase carries an absolute uncertainty
%   of about |U|*eps.  Full precision holds for |U| up to ~1e12; beyond
%   that the error grows linearly and by |U| ~ 1e16 the phase is lost
%   entirely.  This bound is shared by every double-precision
%   implementation (including MATLAB's and SciPy's own ELLIPJ).
%
%   Some definitions of the Jacobi elliptic functions use the modulus
%   k instead of the parameter m.  They are related by m = k^2.
%
%   See also ELLIPKE.

%   L. Shure 1-9-88
%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 5.14 $  $Date: 2001/04/15 12:01:40 $
%
%   Modified by Moiseev Igor,
%   moiseev.igor[at]gmail.com
%   34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
%   Date: 2005/10/04
%
%   The modification of orginal script fixes the problem of slow convergence
%   of the AGM algorithm for the value of parameter M=1.


%   ELLIPJ uses the method of the arithmetic-geometric mean
%   described in [1].
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 16-17.6.


if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) || ~isreal(m)
    error('Input arguments must be real. Use ELLIPJI for complex arguments')
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

% GPU dispatch: move computation to GPU if enabled and available
N_el = numel(u);
if has_gpu()
    [sn,cn,dn,am] = gpu_ellipj(u, m, tol);
    return;
end

% Parallel dispatch: split across workers for large inputs
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N_el >= minChunk
    [sn,cn,dn,am] = parallel_ellipj(u, m, tol, nWorkers, minChunk);
    return;
end

am = zeros(size(u));
cn = zeros(size(u));
sn = cn;
dn = sn;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) || any(m > 1),
  error('M must be in the range 0 <= M <= 1.');
end

I = find(m ~= 1 & m ~= 0);
if ~isempty(I)
    % Preserve distinct parameters exactly; tolerance grouping is a silent
    % data substitution and is particularly damaging near m=1.
    m_vals = m(I);
    [mu, ~, K] = unique(m_vals);
    K = K(:).';

    mumax = length(mu);

    % pre-allocate space and augment if needed
	chunk = 7;
	a = zeros(chunk,mumax);
	c = a;
	b = a;
	a(1,:) = ones(1,mumax);
	c(1,:) = sqrt(mu);
	b(1,:) = sqrt(1-mu);
	n = zeros(1,mumax);
	i = 1;
	while any(abs(c(i,:)) > tol)                                    % Arithmetic-Geometric Mean of A, B and C
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,mumax)];
          b = [b; zeros(2,mumax)];
          c = [c; zeros(2,mumax)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        mask = (abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol);
        n(mask) = i-1;
	end

    mmax = length(I);
    % Use the platform complete integral for reduction.  Deriving K from the
    % tolerance-stopped AGM is adequate for small u but its last-bit error is
    % multiplied by the period count for large u, especially near m=1.
    K_unique = carlsonRF(zeros(size(mu)), 1-mu, ones(size(mu)));
    K_vals = K_unique(K);
    period = floor((u(I) + K_vals) ./ (2 .* K_vals));
    u_reduced = u(I) - 2 .* period .* K_vals;
	phin = zeros(1,mmax);
	a_fin = reshape(a(sub2ind(size(a), n(K) + 1, K)), 1, mmax);   % per-element converged AGM row (K may be a column)
	phin(:) = (2 .^ n(K)).*a_fin.*u_reduced;
	while i > 1
        i = i - 1;
        mask = n(K) >= i;
        if any(mask)
          % asin(c sin/a) = atan2(c sin, sqrt(a^2 cos^2 + b^2 sin^2)) using
          % a^2 - c^2 = b^2: no asin near +/-1, which lost ~7 digits as m -> 1.
          sp = sin(phin(mask));  cp = cos(phin(mask));
          phin(mask) = 0.5*(atan2(c(i+1,K(mask)).*sp, ...
                sqrt((a(i+1,K(mask)).*cp).^2 + (b(i+1,K(mask)).*sp).^2)) + phin(mask));
        end
	end
    quasi_sign = 1 - 2 .* mod(period, 2);
	am(I) = phin + period .* pi;
	sn(I) = quasi_sign .* sin(phin);
	cn(I) = quasi_sign .* cos(phin);
	dn(I) = sqrt((1 - m(I)) + m(I).*cn(I).^2);
end

% Special cases: m = {0, 1}
m0 = find(m == 0);
am(m0) = u(m0);
sn(m0) = sin(u(m0));
cn(m0) = cos(u(m0));
dn(m0) = ones(size(m0));

m1 = find(m == 1);
am(m1) = asin(tanh(u(m1)));
sn(m1) = tanh(u(m1));
cn(m1) = sech(u(m1));
dn(m1) = sech(u(m1));


function [sn,cn,dn,am] = parallel_ellipj(u, m, tol, nWorkers, minChunk)
%PARALLEL_ELLIPJ  Internal helper: split work across parfor workers.
    origSize = size(u);
    N = numel(u);
    u_flat = u(:).'; m_flat = m(:).';
    nChunks = min(nWorkers, ceil(N / minChunk));
    chunkSize = ceil(N / nChunks);
    sn_c = cell(1, nChunks); cn_c = cell(1, nChunks);
    dn_c = cell(1, nChunks); am_c = cell(1, nChunks);
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
            repmat({'ellipj'}, 1, nChunks), u_chunks, m_chunks, tol_c, ...
            'UniformOutput', false);
        for w = 1:nChunks
            sn_c{w} = results{w}{1}; cn_c{w} = results{w}{2};
            dn_c{w} = results{w}{3}; am_c{w} = results{w}{4};
        end
    else
        parfor w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            [sn_c{w}, cn_c{w}, dn_c{w}, am_c{w}] = ellipj(u_flat(i1:i2), m_flat(i1:i2), tol);
        end
    end
    sn = reshape([sn_c{:}], origSize);
    cn = reshape([cn_c{:}], origSize);
    dn = reshape([dn_c{:}], origSize);
    am = reshape([am_c{:}], origSize);


function [sn,cn,dn,am] = gpu_ellipj(u, m, tol)
%GPU_ELLIPJ  Internal helper: compute ellipj using gpuArray.
%   Compatible with both MATLAB gpuArray and Octave ocl package.
    origSize = size(u);
    am = zeros(origSize);
    cn = zeros(origSize);
    sn = zeros(origSize);
    dn = zeros(origSize);

    m = m(:);
    u = u(:);

    if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

    I = find(m ~= 1 & m ~= 0);
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

        % Reduce by the 2K quasi-period before the amplified Landen phase.
        % This mirrors the serial path and prevents large-argument phase loss.
        a_cpu = gather(a);  a_final = a_cpu(sub2ind(size(a_cpu), (1:mmax)', n + 1));   % per-element converged row
        K_vals = carlsonRF(zeros(size(mu)), 1-mu, ones(size(mu)));
        period = floor((u(I) + K_vals) ./ (2 .* K_vals));
        u_reduced = u(I) - 2 .* period .* K_vals;

        % Ascending Landen back-substitution with multiplicative masking
        phin = gpuArray((2 .^ n) .* a_final .* u_reduced);
        for jj = ii-1:-1:1
            active = gpuArray(double(n >= jj));
            sp = sin(phin);  cp = cos(phin);
            phin_new = 0.5*(atan2(c(:,jj+1).*sp, sqrt((a(:,jj+1).*cp).^2 + (b(:,jj+1).*sp).^2)) + phin);
            phin = phin + active .* (phin_new - phin);
        end

        quasi_sign = 1 - 2 .* mod(period, 2);
        phin_cpu = gather(phin);
        cn_val = quasi_sign .* cos(phin_cpu);      % keep everything a column:
        am(I) = phin_cpu + period .* pi;            % m(I) is a column here but
        sn(I) = quasi_sign .* sin(phin_cpu);        % cn(I) indexes a row, so
        cn(I) = cn_val;                             % m(I).*cn(I) was an outer
        dn(I) = sqrt((1 - m(I)) + m(I) .* cn_val.^2);   % product (found on L4)
    end

    % Special cases: m = {0, 1}
    m0 = find(m == 0);
    am(m0) = u(m0);  sn(m0) = sin(u(m0));  cn(m0) = cos(u(m0));  dn(m0) = 1;

    m1 = find(m == 1);
    am(m1) = asin(tanh(u(m1)));  sn(m1) = tanh(u(m1));
    cn(m1) = sech(u(m1));        dn(m1) = sech(u(m1));
