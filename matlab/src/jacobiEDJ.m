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
if compute_J
    [~, ~, ~, phi] = ellipj(u, m);
    [B, D, J] = ellipticBDJ(phi, m, n);
    Eu = u - m .* D;    % E_u = u - m*D_u
    Du = D;
    Ju = J;
else
    [~, ~, ~, phi] = ellipj(u, m);
    [~, D] = ellipticBDJ(phi, m);
    Eu = u - m .* D;
    Du = D;
    Ju = [];
end
