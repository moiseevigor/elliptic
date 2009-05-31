function Pi = elliptic3(u,m,c);
% ELLIPTIC3 evaluates incomplete elliptic integral of the third kind.
%   Pi = ELLIPTIC3(U,M,C) where U is a phase in radians, 0<M<1 is 
%   the module and 0<C<1 is a parameter. 
%
%   ELLIPTIC3 uses Gauss-Legendre 10 points quadrature template 
%   described in [3] to determine the value of the Incomplete Elliptic 
%   Integral of the Third Kind (see [1, 2]).
%
%   Pi(u,m,c) = int(1/((1 - c*sin(t)^2)*sqrt(1 - m*sin(t)^2)), t=0..u)
%
%   Tables generating code ([1], pp. 625-626):
%	    [phi,alpha,c] = meshgrid(0:15:90, 0:15:90, 0:0.1:1); 
%   	Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
%  
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 17.7.
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989.
%   [3] S. Zhang, J. Jin "Computation of Special Functions" (Wiley, 1996).

%   For support, please reply to 
%       moiseev[at]sissa.it
%       Moiseev Igor, 
%       34106, SISSA, via Beirut n. 2-4,  Trieste, Italy

if nargin<3, error('Not enough input arguments.'); end
if ~isreal(u) | ~isreal(m) | ~isreal(c)
    error('Input arguments must be real.')
end
if any(m < 0) | any(m > 1) | any(c < 0) | any(c > 1),  
  error('M and C must be in the range [0, 1].');
end
if any(u > pi/2) | any(u < 0),  
    error('U must be in the range [0, pi/2].'); 
end

[mm,nm] = size(m);
[mu,nu] = size(u);
if length(m)==1, m = m(ones(size(u))); end
if length(c)==1, c = c(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m), size(c), size(u)), 
        error('U, M and C must be the same size.'); 
end

Pi = zeros(size(u));
m = m(:).';    % make a row vector
u = u(:).';
c = c(:).';

I = find( u==pi/2 & m==1 | u==pi/2 & c==1 );

t = [ 0.9931285991850949,  0.9639719272779138,...            % Base points 
      0.9122344282513259,  0.8391169718222188,...            % for Gauss-Legendre integration
      0.7463319064601508,  0.6360536807265150,...
      0.5108670019508271,  0.3737060887154195,...
      0.2277858511416451,  0.07652652113349734 ];                             
w = [ 0.01761400713915212, 0.04060142980038694,...           % Weights
      0.06267204833410907, 0.08327674157670475,...           % for Gauss-Legendre integration
      0.1019301198172404,  0.1181945319615184,...
      0.1316886384491766,  0.1420961093183820,...
      0.1491729864726037,  0.1527533871307258  ];
  
P = 0;  i = 0;
while i < 10
    i  = i + 1;
    c0 = u.*t(i)/2;
    P  = P + w(i).*(g(u/2+c0,m,c) + g(u/2-c0,m,c));
end
P = u/2.*P;
Pi(:) = P;                                                   % Incomplete elliptic integral of the third kind

% special values u==pi/2 & m==1 | u==pi/2 & c==1
Pi(I) = inf;
return;


function g = g(u,m,c)
%  g = 1/((1 - c*sin(u)^2)*sqrt(1 - m*sin(u)^2));

 sn2 = sin(u).^2;
 g = 1./((1 - c.*sn2).*sqrt(1 - m.*sn2));
return;