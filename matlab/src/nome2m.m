function m = nome2m(q)
%NOME2M  Inverse of NOMEQ:  q -> m  (0 < q < 1).  Alias of INVERSENOMEQ.
%   M = NOME2M(Q) returns the parameter m whose nome is Q, elementwise, via
%   the closed theta form m = (theta2(0,q)/theta3(0,q))^4 (DLMF 20.9.1), which
%   is exact on the whole open interval.  The previous fzero bracket
%   [1e-8, 1-1e-8] covered only q in (6e-10, 0.62), and its objective captured
%   the whole input array, so any array input errored inside fzero.
%
%   See also INVERSENOMEQ, NOMEQ.
if ~isreal(q) || any(~(q > 0 & q < 1))
    error('nome2m: q must satisfy 0 < q < 1.');
end
m = inversenomeq(q);
end
