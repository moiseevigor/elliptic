function m = inversenomeq(q)
%INVERSENOMEQ gives the value of Nome m = m(q).
%
%   M = inversenomeq(q), where Q is the Nome of q-series.
%
%   WARNING. The function INVERSENOMEQ does not return correct
%   values of M for Q > 0.6, because of computer precision limitation.
%   The function NomeQ(m) has an essential singularity at M = 1, so
%   it cannot be inverted at this point and actually is very hard to
%   fing and inverse in the neigborhood also.
%   More preciesly:
%        nomeq(1) = 1
%        nomeq(1-eps) = 0.77548641878026
%
%   Example:
%       nomeq(inversenomeq([0.001 0.3 0.4 0.5 0.6 0.7 0.8]))
%
%   See also
%        Standard: ELLIPKE, ELLIPJ
%        Moiseev's package: ELLIPTIC12I, ELLIPTIC3, THETA, AGM.
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


if nargin<1, error('Not enough input arguments.'); end

if ~isreal(q)
    error('Input arguments must be real.')
end

m = zeros(size(q));
q = q(:).';    % make a row vector

bad = isnan(q);
if ~all(q(~bad) >= 0) || ~all(q(~bad) < 1)
    error('Input arguments must be from the interval [0,1).')
end
q(bad) = 0;                         % computed as m(0) = 0, overwritten with NaN below

% Closed form, DLMF 20.9.1:  m = (theta2(0,q)/theta3(0,q))^4
%   theta2(0,q) = 2*q^(1/4) * sum q^(n(n+1)),  theta3(0,q) = 1 + 2*sum q^(n^2)
% The q^(1/4) factor is kept outside the ratio so tiny q cannot underflow.
% This replaces the old interpolation tables, which were documented as
% unreliable for q < 1e-5 and q > 0.76; the series is exact at every scale
% down to m(1e-30) = 1.6e-29.
% Above q_max = 0.778953424877990 the true 1-m = m(exp(pi^2/ln q)) ~ 16 exp(-pi^2/ln(1/q))
% is below eps/2, so the correctly rounded double is exactly 1; the 30-term
% series is not converged there and returned m > 1 (1.034 at q = 0.999).
q_max = 0.778953424877990;
s2 = ones(size(q));                 % sum q^(n(n+1)), n >= 0
s3 = ones(size(q));                 % theta3 = 1 + 2*sum q^(n^2)
qs = min(q, q_max);
for n = 1:30
    s2 = s2 + qs.^(n*(n+1));
    s3 = s3 + 2*qs.^(n^2);
end
m(:) = min(16*qs .* (s2./s3).^4, 1);
m(q > q_max) = 1;
m(bad) = NaN;

% END FUNCTION inversenomeq()
