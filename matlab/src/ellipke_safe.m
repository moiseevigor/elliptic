function [K, E] = ellipke_safe(m, tol)
%ELLIPKE_SAFE  Complete elliptic integrals K(m), E(m) with NaN propagation.
%   [K, E] = ELLIPKE_SAFE(M) is ELLIPKE(M) except that NaN elements of M
%   give NaN instead of aborting: Octave's ellipke raises "algorithm did not
%   converge" as soon as one element is NaN, which took every theta and nome
%   function down with it.  TOL is passed through when given.
if nargin < 2, tol = eps; end
K = nan(size(m));  E = K;
ok = ~isnan(m);
if any(ok(:))
    [K(ok), E(ok)] = ellipke(m(ok), tol);
end
end
