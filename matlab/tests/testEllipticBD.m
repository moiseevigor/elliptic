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

% ------------------------------------------------------------------
%  Group 2: Hardcoded reference table
%
%  Reference values computed by evaluating B=(E-(1-m)*K)/m, D=(K-E)/m,
%  S=(D-B)/m using K,E from ellipke (Gauss AGM, machine epsilon).
%  Cross-checked against carlsonRD via RD(0,1-m,1)/3 = D(m).
%  m = 0 limit is exact: B=D=pi/4, S=pi/16.
% ------------------------------------------------------------------

%!test
%! % ellipticBD reference table at m = 0.0..0.9
%! m_v = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%! expectedB = [ ...
%!   0.7853981633974483, 0.7956042304956568, 0.8066808960371519, ...
%!   0.8188015022917049, 0.8322012900067004, 0.8472130847939785, ...
%!   0.8643348918741709, 0.8843737533686884, 0.9088110737045851, ...
%!   0.9410728015213956];
%! expectedD = [ ...
%!   0.7853981633974483, 0.8168371182245626, 0.8529427025733760, ...
%!   0.8950879458870861, 0.9453180814845530, 1.006861592507393,  ...
%!   1.085232857931855,  1.190989381923781,  1.348394253116269,  ...
%!   1.637019311826777];
%! expectedS = [ ...
%!   0.1963495408493621, 0.2123288772890586, 0.2313090326811207, ...
%!   0.2542881453179373, 0.2827919786946315, 0.3192970154268293, ...
%!   0.3681632767628066, 0.4380223265072747, 0.5494789742646045, ...
%!   0.7732739003393130];
%! [B, D, S] = ellipticBD(m_v);
%! assert(norm(B - expectedB) < 1e-13, 'B reference table failed');
%! assert(norm(D - expectedD) < 1e-13, 'D reference table failed');
%! assert(norm(S - expectedS) < 1e-13, 'S reference table failed');
