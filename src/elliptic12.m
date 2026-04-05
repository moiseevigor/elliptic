function [F,E,Z] = elliptic12(u,m,tol)
% ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals
% of the First, Second Kind and Jacobi's Zeta Function.
%
%   [F,E,Z] = ELLIPTIC12(U,M,TOL) where U is a phase in radians, 0<M<1 is
%   the module and TOL is the tolerance (optional). Default value for
%   the tolerance is eps = 2.220e-16.
%
%   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean
%   and Descending Landen Transformation described in [1] Ch. 17.6,
%   to determine the value of the Incomplete Elliptic Integrals
%   of the First, Second Kind and Jacobi's Zeta Function [1], [2].
%
%       F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
%       E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
%       Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
%
%   Tables generating code ([1], pp. 613-621):
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees
%       [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC12I, ELLIPTIC3, THETA, AGM.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989

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
%
% The code is optimized for ordered inputs produced by the functions
% meshgrid, ndgrid. To obtain maximum performace (up to 30%) for singleton,
% 1-dimensional and random arrays remark call of the function unique(.)
% and edit further code.

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) || ~isreal(m)
    error('Input arguments must be real. Use ELLIPTIC12i for complex arguments.');
end

if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

% Parallel dispatch: split across workers for large inputs
N = numel(u);
nWorkers = get_nworkers();
minChunk = elliptic_config('chunk_size');
if nWorkers > 1 && N >= minChunk
    [F,E,Z] = parallel_elliptic12(u, m, tol, nWorkers, minChunk);
    return;
end

F = zeros(size(u));
E = F;
Z = E;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

% check whether we've been asked to evaluate the integrals for values
% smaller than eps = 2.220446049250313e-16, if so we suppose it equal zero
m(m<eps) = 0;

I = uint32( find(m ~= 1 & m ~= 0) );
if ~isempty(I)
    % Use standard uniquetol for numerical precision issues
    % This is the recommended MATLAB approach since R2015a
    m_vals = m(I);
    tol_unique = 1e-11;

    [mu, ~, K] = uniquetol_compat(m_vals, tol_unique);
    K = uint32(K(:).');  % Ensure K is a row vector
    mumax = length(mu);
    signU = sign(u(I));

    % pre-allocate space and augment if needed
	chunk = 7;
	a = zeros(chunk,mumax);
	c = a;
	b = a;
	a(1,:) = ones(1,mumax);
	c(1,:) = sqrt(mu);
	b(1,:) = sqrt(1-mu);
	n = uint32( zeros(1,mumax) );
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
	mn = double(max(n));
	phin = zeros(1,mmax);     C  = zeros(1,mmax);
	Cp = C;  e  = zeros(1,mmax);  phin(:) = signU.*u(I);
	c2 = c.^2;
	e_vals = 2.^(0:mn-1);                                          % pre-compute powers of 2
	for i = 1:mn                                                    % Descending Landen Transformation
        mask = n(K) > i;
        if any(mask)
            phin(mask) = atan(b(i,K(mask))./a(i,K(mask)).*tan(phin(mask))) + ...
                pi.*ceil(phin(mask)/pi - 0.5) + phin(mask);
            e(mask) = e_vals(i);
            C(mask) = C(mask)  + double(e_vals(i))*c2(i,K(mask));
            Cp(mask)= Cp(mask) + c(i+1,K(mask)).*sin(phin(mask));
        end
	end

    Ff = phin ./ (a(mn,K).*double(e)*2);
    F(I) = Ff.*signU;                                               % Incomplete Ell. Int. of the First Kind
    Z(I) = Cp.*signU;                                               % Jacobi Zeta Function
    E(I) = (Cp + (1 - 1/2*C) .* Ff).*signU;                         % Incomplete Ell. Int. of the Second Kind
end

% Special cases: m == {0, 1}
m0 = find(m == 0);
if ~isempty(m0), F(m0) = u(m0); E(m0) = u(m0); Z(m0) = 0; end

m1 = find(m == 1);
um1 = abs(u(m1));
if ~isempty(m1),
    N = floor( (um1+pi/2)/pi );
    M = find(um1 < pi/2);

    F(m1(M)) = log(tan(pi/4 + u(m1(M))/2));
    F(m1(um1 >= pi/2)) = Inf.*sign(u(m1(um1 >= pi/2)));

    E(m1) = ((-1).^N .* sin(um1) + 2*N).*sign(u(m1));

    Z(m1) = (-1).^N .* sin(u(m1));
end


function [F,E,Z] = parallel_elliptic12(u, m, tol, nWorkers, minChunk)
%PARALLEL_ELLIPTIC12  Internal helper: split work across parfor workers.
    origSize = size(u);
    N = numel(u);
    u_flat = u(:).'; m_flat = m(:).';
    nChunks = min(nWorkers, ceil(N / minChunk));
    chunkSize = ceil(N / nChunks);
    F_cells = cell(1, nChunks);
    E_cells = cell(1, nChunks);
    Z_cells = cell(1, nChunks);
    if exist('OCTAVE_VERSION', 'builtin')
        u_chunks = cell(1, nChunks);
        m_chunks = cell(1, nChunks);
        for w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            u_chunks{w} = u_flat(i1:i2);
            m_chunks{w} = m_flat(i1:i2);
        end
        tol_c = repmat({tol}, 1, nChunks);
        results = parcellfun(nWorkers, @par_worker, ...
            repmat({'elliptic12'}, 1, nChunks), u_chunks, m_chunks, tol_c, ...
            'UniformOutput', false);
        for w = 1:nChunks
            F_cells{w} = results{w}{1};
            E_cells{w} = results{w}{2};
            Z_cells{w} = results{w}{3};
        end
    else
        parfor w = 1:nChunks
            i1 = (w-1)*chunkSize + 1;
            i2 = min(w*chunkSize, N);
            [F_cells{w}, E_cells{w}, Z_cells{w}] = elliptic12(u_flat(i1:i2), m_flat(i1:i2), tol);
        end
    end
    F = reshape([F_cells{:}], origSize);
    E = reshape([E_cells{:}], origSize);
    Z = reshape([Z_cells{:}], origSize);
