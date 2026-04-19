% Tests for Bulirsch's generalised complete elliptic integral cel(kc,p,a,b)
% and thin wrappers cel1, cel2, cel3.
%
% Run with:  test testBulirschCEL

% ------------------------------------------------------------------
%  Group 1: Special cases matching K, E, B, D, Pi
% ------------------------------------------------------------------

%!test
%! % cel1(kc) = K(1-kc^2)
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! kc = sqrt(1 - m);
%! [K, ~] = ellipke(m);
%! C1 = cel1(kc);
%! assert(max(abs(C1 - K)) < 1e-13, 'cel1(kc) != K(m)');

%!test
%! % cel2(kc, 1, kc^2) = E(m)  [since kc^2 = 1-m]
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! kc = sqrt(1 - m);
%! [~, E] = ellipke(m);
%! C2 = cel2(kc, 1, kc.^2);    % cel(kc, 1, 1, kc^2) = a*B + b*D = B + kc^2*D = E
%! assert(max(abs(C2 - E)) < 1e-13, 'cel2(kc,1,kc^2) != E(m)');

%!test
%! % cel2(kc, 1, 0) = B(m)
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! kc = sqrt(1 - m);
%! [B, ~, ~] = ellipticBD(m);
%! C2B = cel2(kc, 1, 0);
%! assert(max(abs(C2B - B)) < 1e-13, 'cel2(kc,1,0) != B(m)');

%!test
%! % cel2(kc, 0, 1) = D(m)
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! kc = sqrt(1 - m);
%! [~, D, ~] = ellipticBD(m);
%! C2D = cel2(kc, 0, 1);
%! assert(max(abs(C2D - D)) < 1e-13, 'cel2(kc,0,1) != D(m)');

%!test
%! % cel3(kc, 1-n) = Pi(n|m)  [complete 3rd kind]
%! m = 0.5;  n = 0.3;
%! kc = sqrt(1 - m);
%! Pi_leg = elliptic3(pi/2, m, n);
%! Pi_cel = cel3(kc, 1-n);
%! assert(abs(Pi_cel - Pi_leg) < 1e-12, ...
%!     sprintf('cel3(kc,1-n) != Pi: %g vs %g', Pi_cel, Pi_leg));

% ------------------------------------------------------------------
%  Group 2: General cel formula versus numerical integration
% ------------------------------------------------------------------

%!test
%! % General p≠1: cel vs direct numerical integration
%! m = 0.7;  kc = sqrt(1-m);
%! phi = linspace(0, pi/2, 100000);
%! cos2 = cos(phi).^2; sin2 = sin(phi).^2;
%! Delta = sqrt(cos2 + kc^2.*sin2);
%! tests = [0.4, 2.0, 0.5;  % [p, a, b]
%!          1.8, 0.7, 1.3;
%!          0.1, 0.3, 1.5;
%!          2.5, 1.0, 0.0];
%! for k = 1:size(tests,1)
%!     p_t = tests(k,1); a_t = tests(k,2); b_t = tests(k,3);
%!     num  = a_t.*cos2 + b_t.*sin2;
%!     denom = (cos2 + p_t.*sin2) .* Delta;
%!     cel_num = trapz(phi, num./denom);
%!     cel_my  = cel(kc, p_t, a_t, b_t);
%!     assert(abs(cel_my - cel_num) < 1e-10, ...
%!         sprintf('cel(kc,%g,%g,%g): diff=%g', p_t, a_t, b_t, abs(cel_my-cel_num)));
%! end

% ------------------------------------------------------------------
%  Group 3: Array inputs and edge cases
% ------------------------------------------------------------------

%!test
%! % Array over kc (different m values)
%! m = linspace(0.1, 0.9, 10);
%! kc = sqrt(1 - m);
%! [K, ~] = ellipke(m);
%! C = cel1(kc);
%! assert(max(abs(C - K)) < 1e-13, 'cel1 array failed');

%!test
%! % cel(kc, p, a, b) is linear in a, b
%! kc = 0.6; p = 0.7; a1 = 1; b1 = 2; a2 = 3; b2 = 0.5;
%! c1 = cel(kc, p, a1, b1);
%! c2 = cel(kc, p, a2, b2);
%! c3 = cel(kc, p, a1+a2, b1+b2);
%! assert(abs(c1 + c2 - c3) < 1e-13, 'cel not linear in a,b');

%!test
%! % kc=0 case: K→∞, cel should diverge (pi/2 * [a*K + ...]→ Inf)
%! % For kc=0 and p=1: cel = a*B(1) + b*D(1), D→∞
%! % Just verify cel is finite for kc≥eps
%! kc = 1e-8; m_near1 = 1 - kc^2;
%! [K, ~] = ellipke(m_near1);
%! C = cel1(kc);
%! assert(isfinite(C), 'cel1 near kc=0 should be finite (large)');
%! assert(C > 0, 'cel1 should be positive');

% ------------------------------------------------------------------
%  Group 4: Hardcoded reference table for cel(kc, p, a, b)
%
%  Values computed via Carlson RJ decomposition:
%    p=1 case:  cel = a*B(m) + b*D(m)
%    p≠1 case:  cel = a*K + (b-a*p)*(Pi(1-p|m)-K)/(1-p)
%  Cross-checked against numerical integration (trapz over [0,pi/2])
%  with agreement to < 1e-10.
%
%  cel1(kc) = K(1-kc^2) — exact to machine epsilon.
%  cel2(kc,1,kc^2) = E(1-kc^2) — exact to machine epsilon.
% ------------------------------------------------------------------

%!test
%! % cel generic reference table: 10 (kc, p, a, b) test points
%! tests = [ ...
%!   0.5, 0.3, 1.0, 0.5,  2.785372539236548; ...
%!   0.5, 1.0, 0.8, 1.2,  2.229457648629679; ...
%!   0.5, 2.0, 0.5, 1.5,  1.436498488170953; ...
%!   0.7, 0.4, 1.0, 1.0,  3.066196685396346; ...
%!   0.7, 1.0, 1.0, 1.0,  1.862640802332739; ...
%!   0.7, 1.5, 0.3, 0.7,  0.7431229068315112; ...
%!   0.3, 0.2, 0.5, 1.0,  6.654375216945862; ...
%!   0.3, 1.0, 1.0, 0.09, 1.096477517392227; ...
%!   0.9, 0.1, 1.0, 1.0,  5.378472453168364; ...
%!   0.9, 1.0, 1.0, 1.0,  1.654616667522527];
%! for k = 1:size(tests, 1)
%!   kc = tests(k,1); p = tests(k,2); a = tests(k,3);
%!   b  = tests(k,4); ref = tests(k,5);
%!   cv = cel(kc, p, a, b);
%!   assert(abs(cv - ref) < 1e-13, ...
%!       sprintf('cel(%.2g,%.2g,%.2g,%.2g): got %.16g, expected %.16g', kc,p,a,b,cv,ref));
%! end

%!test
%! % cel1 reference table: kc in {0.2, 0.4, 0.6, 0.7071, 0.8944}
%! kc_v = [0.2, 0.4, 0.6, 1/sqrt(2), sqrt(1-0.1)];
%! m_v  = 1 - kc_v.^2;
%! [Kv, ~] = ellipke(m_v);
%! C1 = cel1(kc_v);
%! assert(max(abs(C1 - Kv)) < 1e-13, 'cel1 reference table failed');

%!test
%! % cel2(kc, 1, kc^2) reference table: equals E(1-kc^2)
%! kc_v = [0.3, 0.5, 0.7, 0.8, 0.9];
%! m_v  = 1 - kc_v.^2;
%! [~, Ev] = ellipke(m_v);
%! C2 = cel2(kc_v, 1, kc_v.^2);
%! assert(max(abs(C2 - Ev)) < 1e-13, 'cel2(kc,1,kc^2)=E reference table failed');
