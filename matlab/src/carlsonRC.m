function RC = carlsonRC(x, y)
%CARLSONRC  Carlson's degenerate symmetric elliptic integral R_C(x,y).
%
%   RC = CARLSONRC(X, Y) evaluates
%
%       R_C(x,y) = (1/2) * ∫_0^∞ (t+x)^(-1/2) (t+y)^(-1) dt
%
%   which equals R_F(x,y,y) (DLMF 19.16.6).
%
%   Closed-form evaluation (DLMF 19.2.17–19.2.19):
%
%       x >= 0, y > 0:
%           y > x : R_C = arctan(√((y-x)/x))  / √(y-x)   (DLMF 19.2.17)
%           y < x : R_C = arctanh(√((x-y)/x)) / √(x-y)   (DLMF 19.2.18)
%           y = x : R_C = 1/√x
%       x = 0    : R_C = π / (2√y)
%
%   Degenerate cases returned as Inf when y=0.
%
%   X and Y may be any combination of scalars or same-size arrays.
%   Scalar inputs are broadcast to the larger size.
%
%   No separate parallel/GPU dispatch is provided at this level; the
%   function is called element-wise from higher-level routines that own
%   the dispatch logic.
%
%   References:
%   [1] NIST DLMF §19.2, §19.16  https://dlmf.nist.gov/19.2
%   [2] B.C. Carlson, "Computing Elliptic Integrals by Duplication,"
%       Numer. Math. 33 (1979), 1–16.

% Empty input -> empty output of the same shape (elementwise semantics; the
% size checks below would otherwise reject [] against a scalar).
if nargin >= 2 && (isempty(x) || isempty(y))
    sz = size(x);
    if isempty(y), sz = size(y); end
    RC = zeros(sz);
    return;
end

if nargin < 2, error('carlsonRC: requires two arguments (x, y).'); end
if ~isreal(x) || ~isreal(y)
    error('carlsonRC: all input arguments must be real.');
end

[x, y] = carlsonRC_broadcast(x, y);
origSize = size(x);
x = x(:).';  y = y(:).';

% NaN, Inf and negative x give NaN (y < 0 is the principal value, handled by the core).
bad = ~(x >= 0) | isinf(x) | isnan(y) | isinf(y);
x(bad) = 1;  y(bad) = 1;
RC = carlsonRC_core(x, y);
RC(bad) = NaN;
RC = reshape(RC, origSize);


% -----------------------------------------------------------------------
function RC = carlsonRC_core(x, y)
%CARLSONRC_CORE  Vectorised serial evaluation (row-vector inputs).
N  = numel(x);
RC = zeros(1, N);

% Degenerate: y = 0  → Inf
pole = (y == 0);

% x = 0: R_C(0,y) = π/(2√y)
zero_x = (x == 0) & ~pole;

% y = x (and x > 0): R_C = 1/√x
eq = (x == y) & ~pole & ~zero_x;

% y > x (and x > 0): arctan branch
gt = (y > x) & ~pole & ~zero_x & ~eq;

% y < x: arctanh branch
lt = (y < x) & ~pole & ~zero_x & ~eq;

RC(pole)  = Inf;
RC(zero_x) = pi ./ (2 .* sqrt(y(zero_x)));
RC(eq)     = 1 ./ sqrt(x(eq));

if any(gt)
    d      = sqrt((y(gt) - x(gt)) ./ x(gt));
    RC(gt) = atan(d) ./ (sqrt(y(gt) - x(gt)));
end

if any(lt)
    % log((sqrt(x)+sqrt(x-y))/sqrt(y))/sqrt(x-y) == atanh(sqrt(1-y/x))/sqrt(x-y)
    % without the 1 - sqrt(1-eps) cancellation (RC(3,1e-10) lost 8 digits).
    % ... as log1p: log((sx+sxy)/sy) = log1p(((x-y)/(sx+sy) + sxy)/sy); the
    % plain log lost 9 digits again for tiny x - y (RC(1+1e-13, 1)).
    xl = x(lt);  yl = y(lt);  sx = sqrt(xl);  sy = sqrt(yl);  sxy = sqrt(xl - yl);
    RC(lt) = log1p(((xl - yl)./(sx + sy) + sxy) ./ sy) ./ sxy;
end


% -----------------------------------------------------------------------
function [x, y] = carlsonRC_broadcast(x, y)
refSz = [1 1];
for v = {x, y}
    if numel(v{1}) > 1, refSz = size(v{1}); break; end
end
if isscalar(x), x = x(ones(refSz)); end
if isscalar(y), y = y(ones(refSz)); end
if ~isequal(size(x), size(y))
    error('carlsonRC: x and y must be the same size or scalar.');
end
