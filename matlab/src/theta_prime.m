function [th, thp] = theta_prime(j, z, m, tol)
%THETA_PRIME evaluates theta functions and their derivatives.
%   [TH, THP] = THETA_PRIME(J, Z, M) returns values of the Jacobi theta
%   function TH and its derivative THP with respect to the argument Z.
%   J is the type of theta function (1, 2, 3, or 4), Z is the argument,
%   and M is the parameter (0 <= M <= 1).
%
%   [TH, THP] = THETA_PRIME(J, Z, M, TOL) computes the theta function
%   and its derivative to the accuracy TOL instead of the default TOL = EPS.
%
%   The arrays Z and M must be the same size (or either can be scalar).
%   As currently implemented, M is limited to 0 <= M <= 1.
%
%   The theta functions are defined as:
%       θ₁(z|τ) = 2∑[n=0,∞] (-1)ⁿ q^((n+1/2)²) sin((2n+1)z)
%       θ₂(z|τ) = 2∑[n=0,∞] q^((n+1/2)²) cos((2n+1)z)
%       θ₃(z|τ) = 1 + 2∑[n=1,∞] q^(n²) cos(2nz)
%       θ₄(z|τ) = 1 + 2∑[n=1,∞] (-1)ⁿ q^(n²) cos(2nz)
%
%   where q = exp(iπτ) is the nome and τ is related to parameter M.
%
%   The derivatives are computed using the relation:
%       θ'ⱼ(z,m) = θⱼ(z,m) * (2K/π) * (Z + δⱼ)
%   where K is the complete elliptic integral, Z is the Jacobi zeta
%   function, and δⱼ depends on the theta function type.
%
%   The parameter M is related to the nome Q as Q = exp(-π*K(1-M)/K(M)).
%   Some definitions of the Jacobi elliptic functions use the modulus
%   k instead of the parameter m. They are related by m = k².
%
%   Example:
%       j = 1; z = 0.5; m = 0.8;
%       [th, thp] = theta_prime(j, z, m);
%       % Verify using numerical differentiation
%       dz = 1e-8;
%       th_plus = theta(j, z + dz, m);
%       th_minus = theta(j, z - dz, m);
%       thp_numerical = (th_plus - th_minus) / (2 * dz);
%
%   Note: To reproduce results in Mathematica, use the inversenomeq function
%   to convert from Mathematica's nome parameter to the parameter m:
%       [th_4, thp_4] = theta_prime(4, z, inversenomeq(0.1));
%
%   See also THETA, ELLIPJ, ELLIPTIC12, ELLIPKE, JACOBITHETAETA, INVERSENOMEQ.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 16-17.6.
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989
%   [3] E. T. Whittaker and G. N. Watson, "A Course of Modern Analysis"
%       Cambridge University Press, 4th ed., 1996, Ch. 21.

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

% Empty input -> empty output of the same shape (elementwise semantics; the
% size checks below would otherwise reject [] against a scalar).
if nargin >= 3 && (isempty(j) || isempty(z) || isempty(m))
    sz = size(j);
    if isempty(z), sz = size(z); end
    if isempty(m), sz = size(m); end
    th = zeros(sz);
    thp = zeros(sz);
    return;
end

if nargin < 4, tol = eps; end
if nargin < 3, error('Not enough input arguments.'); end

if ~isreal(z) || ~isreal(m)
    error('Input arguments must be real.');
end

if ~isscalar(j) || ~ismember(j, [1, 2, 3, 4])
    error('J must be a scalar integer: 1, 2, 3, or 4.');
end

if length(m) == 1 && length(z) > 1
    m = m(ones(size(z)));
elseif length(z) == 1 && length(m) > 1
    z = z(ones(size(m)));
end

if ~isequal(size(m), size(z))
    error('Z and M must be the same size.');
end

if any(m < 0) || any(m > 1)
    error('M must be in the range 0 <= M <= 1.');
end

K  = ellipke_safe(m);
Kp = carlsonRF(zeros(size(m)), m, ones(size(m)));   % K(1-m) from the exact m (see NOMEQ)
q  = exp(-pi .* Kp ./ K);
q(~(q < 1)) = 0;                  % m == 1 guard
[th, thp] = theta_series(j, z, q, tol);
th(m == 1 | isnan(m)) = NaN;  thp(m == 1 | isnan(m)) = NaN;   % q = 0 stood in for NaN m
if j == 1, th(m == 0) = 0; end

end
