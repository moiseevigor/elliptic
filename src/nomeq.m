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
%     moiseev[at]sissa.it, moiseev.igor[at]gmail.com
%     Moiseev Igor, 
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy


if nargin<2, tol = eps; end
if nargin<1, error('Not enough input arguments.'); end

if ~isreal(m)
    error('Input arguments must be real.')
end

NomeQ = exp(-pi*ellipke(1-m,tol)./ellipke(m,tol));

% END FUNCTION nomeq()