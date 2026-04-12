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
