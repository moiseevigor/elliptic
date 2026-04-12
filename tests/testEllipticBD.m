% Tests for complete associate elliptic integrals B(m), D(m), S(m).
%
% Run with:  test testEllipticBD

% ------------------------------------------------------------------
%  Group 1: Basic identities relating B, D to K and E
% ------------------------------------------------------------------

%!test
%! % K = B + D
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! [K, E] = ellipke(m);
%! [B, D, ~] = ellipticBD(m);
%! assert(max(abs(B + D - K)) < 1e-14, 'K = B+D failed');

%!test
%! % E = B + mc*D  (mc = 1-m)
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! [K, E] = ellipke(m);
%! [B, D, ~] = ellipticBD(m);
%! mc = 1 - m;
%! assert(max(abs(B + mc.*D - E)) < 1e-14, 'E = B+mc*D failed');

%!test
%! % S = (D-B)/m exactly
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! [~, ~, S] = ellipticBD(m);
%! [~, D, ~] = ellipticBD(m);
%! B_v = ellipticBD(m);
%! assert(max(abs(S - (D - B_v)./m)) < 1e-14, 'S = (D-B)/m failed');

%!test
%! % m=0 limits: B=D=pi/4, S=pi/16
%! [B0, D0, S0] = ellipticBD(0);
%! assert(abs(B0 - pi/4) < 1e-14, 'B(0) != pi/4');
%! assert(abs(D0 - pi/4) < 1e-14, 'D(0) != pi/4');
%! assert(abs(S0 - pi/16) < 1e-14, 'S(0) != pi/16');

%!test
%! % Array output size
%! m = linspace(0, 0.99, 50);
%! [B, D, S] = ellipticBD(m);
%! assert(isequal(size(B), size(m)), 'B size mismatch');
%! assert(isequal(size(D), size(m)), 'D size mismatch');
%! assert(all(isfinite(B)), 'B has non-finite values');
%! assert(all(isfinite(D)), 'D has non-finite values');

%!test
%! % D > 0 and B > 0 for all m in [0,1)
%! m = linspace(0, 0.99, 50);
%! [B, D, ~] = ellipticBD(m);
%! assert(all(B > 0), 'B should be positive');
%! assert(all(D > 0), 'D should be positive');

%!test
%! % Scalar input
%! [B, D, S] = ellipticBD(0.5);
%! [K, E] = ellipke(0.5);
%! assert(abs(B+D-K) < 1e-14, 'K=B+D scalar failed');
%! assert(abs(B+0.5*D-E) < 1e-14, 'E=B+mc*D scalar failed');
