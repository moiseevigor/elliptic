function NomeQ = nomeq(m,tol)
%NOMEQ gives the value of Nome q = q(m).
%
%   NomeQ = nomeq(M,TOL), where 0<=M<=1 is the module and
%   TOL is the tolerance (optional). Default value for
%   the tolerance is eps = 2.220e-16.
%
%   See also
%        Standard: ELLIPKE, ELLIPJ,
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


if nargin<2, tol = eps; end
if nargin<1, error('Not enough input arguments.'); end

if ~isreal(m)
    error('Input arguments must be real.')
end

% K'(m) = K(1-m) = R_F(0, m, 1) evaluated from the EXACT argument m:
% ellipke(1-m) rounds 1-m first and lost ~eps/m relative digits
% (q(1e-16) was 11% off, q(1e-17) came back 0).
% NaN elements propagate; finite elements outside [0, 1] are a domain error
% (ellipke used to abort with 'algorithm did not converge' on a single NaN).
bad = isnan(m);
if any(m(~bad) < 0) || any(m(~bad) > 1)
    error('nomeq: m must be in the range 0 <= m <= 1.');
end
NomeQ = nan(size(m));
mv = m(~bad);
NomeQ(~bad) = exp(-pi*carlsonRF(zeros(size(mv)), mv, ones(size(mv)))./ellipke(mv,tol));

% END FUNCTION nomeq()