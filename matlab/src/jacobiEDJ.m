function [Eu, Du, Ju] = jacobiEDJ(u, m, n)
%JACOBIEDJ  Jacobi-argument associate integrals E_u(u|m), D_u(u|m), J_u(u,n|m).
%
%   [Eu, Du, Ju] = JACOBIEDJ(U, M, N) evaluates the associate integrals
%   expressed in the Jacobi argument u = F(am(u|m) | m):
%
%       E_u(u|m)    = ∫₀^u dn²(v|m) dv  =  E(am(u|m)|m)
%       D_u(u|m)    = ∫₀^u sn²(v|m) dv  =  D(am(u|m)|m)
%       J_u(u,n|m)  = ∫₀^u sn²(v)/( 1−n·sn²(v) ) dv  =  J(am(u|m), n, m)
%
%   [Eu, Du] = JACOBIEDJ(U, M) computes only Eu and Du (Ju = []).
%
%   Relations to the Jacobi argument:
%
%       B_u(u|m) = u − D_u(u|m)   (since B_u + D_u = u)
%       E_u(u|m) = u − m·D_u(u|m) (since E = B + mc·D = (u−D) + mc·D = u − m·D)
%
%   Algorithm: convert Jacobi argument to amplitude φ = am(u|m) via
%   ELLIPJ, then delegate to ELLIPTICBDJ(φ, m, n).
%
%   U, M, N may be scalars or arrays of the same size. M must satisfy
%   0 ≤ m < 1 and |u| ≤ K(m) for real amplitudes (outside this range the
%   functions continue via periodicity of ellipj).
%
%   References:
%   [1] NIST DLMF §22.16  https://dlmf.nist.gov/22.16
%   [2] T. Fukushima, "Elliptic functions and elliptic integrals for
%       celestial mechanics and dynamical astronomy," (2015), §5.

compute_J = (nargin >= 3);

if nargin < 2, error('jacobiEDJ: requires at least two arguments (u, m).'); end
if ~isreal(u) || ~isreal(m)
    error('jacobiEDJ: u and m must be real.');
end
if compute_J && ~isreal(n)
    error('jacobiEDJ: n must be real.');
end

% Get amplitude phi = am(u|m) via ellipj
% ellipj returns [sn, cn, dn, am]
% Reduce u by the period 2K BEFORE taking the amplitude.  am(u) at |u| ~ 1e3
% carries eps*|am| ~ 5e-14 of rounding, and near phi = pi/2 with m -> 1 the
% map phi -> D(phi|m) is steep (dD/dphi = sin^2/Delta ~ 1/sqrt(1-m)), so
% D_u(1520|1-1e-8) was off by 3e-10.  As functions of u the integrals are
% perfectly conditioned (dD_u/du = sn^2 <= 1): reduce u, take the amplitude
% of the reduced argument (|am_r| <= pi/2), add 2k times the complete
% integrals.
if isscalar(m) && ~isscalar(u), m = m + zeros(size(u)); end
if isscalar(u) && ~isscalar(m), u = u + zeros(size(m)); end
if compute_J && isscalar(n) && ~isscalar(u), n = n + zeros(size(u)); end
one = ones(size(m));  zed = zeros(size(m));
K   = carlsonRF(zed, 1 - m, one);
k   = floor((u + K) ./ (2 .* K));
u_r = u - 2 .* k .* K;
[~, ~, ~, phi_r] = ellipj(u_r, m);
D_cpl = carlsonRD(zed, 1 - m, one) ./ 3;               % D(m)
if compute_J
    [~, D_r, J_r] = ellipticBDJ(phi_r, m, n);
    J_cpl = zed;
    kk = find(k ~= 0);
    if ~isempty(kk)
        J_cpl(kk) = carlsonRJ(zed(kk), 1 - m(kk), one(kk), 1 - n(kk)) ./ 3;   % J(n|m)
    end
    Ju = J_r + 2 .* k .* J_cpl;
else
    [~, D_r] = ellipticBDJ(phi_r, m);
    Ju = [];
end
Du = D_r + 2 .* k .* D_cpl;
Eu = u - m .* Du;    % E_u = u - m*D_u  (error eps*|u|, the conditioning floor)
