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
%     moiseev[at]sissa.it, moiseev.igor[at]gmail.com
%     Moiseev Igor, 
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy


if nargin<1, error('Not enough input arguments.'); end

if ~isreal(q)
    error('Input arguments must be real.')
end

m = zeros(size(q));
q = q(:).';    % make a row vector
maxq = max(q);

if ~all(q >= 0) || ~all(q <= 1)
    error('Input arguments must be from the interval [0,1].')
end

if any(q > 0.76) || any(q < 0.00001)
    warning('WarnTests:convertTest', ...
            'The function INVERSENOMEQ does not return \ncorrect values of M for Q < 0.00001 and Q > 0.76, because of computer precision limitation.');
end

I = find (q <= 0.4);
J = find (q > 0.4 & q <= 0.6);
P = find (q > 0.6);

if (~isempty(I))
    mm = 0:0.0001:1;
    K = ellipke(mm);
    KK = K(end:-1:1)./K;
    m(I) = interp1(KK, mm, -1/pi*log(q(I)), 'pchip','extrap');
end

if (~isempty(J))
    mm = 0.9996:0.0000001:1-eps;
    %K = 1/8*(-2+2*mm-2*(-5+mm)*log(4)+(-5+mm).*log(1-mm));
    K = 1/128*(-53+74*mm-21*mm.^2+2*(89+mm.*(-34+9*mm))*log(4)+(-89+mm.*(34-9*mm)).*log(1-mm));
    KK = pi/2./K;
    m(J) = interp1(KK, mm, -1/pi*log(q(J)), 'pchip','extrap');
end

if (~isempty(P))
    mm = (1-10^8*eps):1000*eps:1-eps;
    K = 1/128*(-53+74*mm-21*mm.^2+2*(89+mm.*(-34+9*mm))*log(4)+(-89+mm.*(34-9*mm)).*log(1-mm));
    KK = pi/2./K;
    m(P) = interp1(KK, mm, -1/pi*log(q(P)), 'pchip','extrap');
    % plot(mm,K*log(q(P))+pi*pi/2,'.');
end

% END FUNCTION inversenomeq()
