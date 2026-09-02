function [th, thp] = theta_series(j, v, q, tol)
%THETA_SERIES  Jacobi theta function theta_j(v, q) and d/dv from the q-series.
%   [TH, THP] = THETA_SERIES(J, V, Q, TOL) evaluates (A&S 16.27)
%       theta_1 = 2 sum (-1)^n q^((n+1/2)^2) sin((2n+1)v)
%       theta_2 = 2 sum        q^((n+1/2)^2) cos((2n+1)v)
%       theta_3 = 1 + 2 sum        q^(n^2) cos(2nv)
%       theta_4 = 1 + 2 sum (-1)^n q^(n^2) cos(2nv)
%   and its v-derivative, for arrays V and Q of the same size (0 <= Q < 1).
%
%   sin/cos of the multiples (2n+1)v and 2nv come from the angle-addition
%   recurrence started at sin v, cos v: forming k*v as a double product
%   rounds by eps*|k v|, which cost 1e-12 at v ~ 1e8 and 1e-8 at v ~ 1e11.
%
%   Shared by THETA, THETA_PRIME, JACOBITHETAETA, WEIERSTRASSZETA/SIGMA.

if nargin < 4, tol = eps; end
qmax = max([q(:); 0]);
if qmax > 0
    nT = min(1000, max(1, ceil(sqrt(log(tol) / log(qmax)))));
else
    nT = 1;
end
s1 = sin(v);  c1 = cos(v);
s2 = 2 .* s1 .* c1;  c2 = 1 - 2 .* s1.^2;          % sin 2v, cos 2v
th = zeros(size(v));  thp = th;
switch j
    case 1
        sk = s1;  ck = c1;
        for n = 0:nT
            qq = (-1)^n .* q.^((n+0.5)^2);
            th  = th  + qq .* sk;
            thp = thp + qq .* (2*n+1) .* ck;
            [sk, ck] = deal(sk.*c2 + ck.*s2, ck.*c2 - sk.*s2);
        end
        th = 2*th;  thp = 2*thp;
    case 2
        sk = s1;  ck = c1;
        for n = 0:nT
            qq = q.^((n+0.5)^2);
            th  = th  + qq .* ck;
            thp = thp - qq .* (2*n+1) .* sk;
            [sk, ck] = deal(sk.*c2 + ck.*s2, ck.*c2 - sk.*s2);
        end
        th = 2*th;  thp = 2*thp;
    case 3
        sk = s2;  ck = c2;  th = ones(size(v));
        for n = 1:nT
            qq = q.^(n^2);
            th  = th  + 2 .* qq .* ck;
            thp = thp - 4*n .* qq .* sk;
            [sk, ck] = deal(sk.*c2 + ck.*s2, ck.*c2 - sk.*s2);
        end
    case 4
        sk = s2;  ck = c2;  th = ones(size(v));
        for n = 1:nT
            qq = (-1)^n .* q.^(n^2);
            th  = th  + 2 .* qq .* ck;
            thp = thp - 4*n .* qq .* sk;
            [sk, ck] = deal(sk.*c2 + ck.*s2, ck.*c2 - sk.*s2);
        end
    otherwise
        error('theta_series: J must be 1, 2, 3, or 4.');
end
