% Tests for incomplete associate elliptic integrals B(phi|m), D(phi|m), J(phi,n|m)
% and Jacobi-form integrals Eu, Du, Ju from jacobiEDJ.
%
% Run with:  test testEllipticBDJ

% ------------------------------------------------------------------
%  Group 1: ellipticBDJ — connection to F, E, Pi
% ------------------------------------------------------------------

%!test
%! % F = B + D (versus elliptic12)
%! phi = linspace(0.05, 1.4, 40);  m = 0.5;
%! [B, D] = ellipticBDJ(phi, m);
%! [F, ~] = elliptic12(phi, m .* ones(size(phi)));
%! assert(max(abs(B + D - F)) < 1e-12, 'F = B+D failed');

%!test
%! % E = B + mc*D (versus elliptic12)
%! phi = linspace(0.05, 1.4, 40);  m = 0.5;
%! [B, D] = ellipticBDJ(phi, m);
%! [~, E] = elliptic12(phi, m .* ones(size(phi)));
%! mc = 1 - m;
%! assert(max(abs(B + mc.*D - E)) < 1e-12, 'E = B+mc*D failed');

%!test
%! % Pi = B + D + n*J (versus elliptic3)
%! phi = 0.8;  m = 0.5;  n = 0.3;
%! [B, D, J] = ellipticBDJ(phi, m, n);
%! Pi = elliptic3(phi, m, n);
%! assert(abs(B + D + n*J - Pi) < 1e-11, ...
%!     sprintf('Pi = B+D+nJ failed: diff=%g', abs(B+D+n*J-Pi)));

%!test
%! % phi=0: B = D = J = 0
%! [B0, D0, J0] = ellipticBDJ(0, 0.5, 0.3);
%! assert(B0 == 0 && D0 == 0 && J0 == 0, 'B=D=J=0 at phi=0 failed');

%!test
%! % phi = pi/2: B and D match ellipticBD
%! m = 0.6;
%! [B_inc, D_inc] = ellipticBDJ(pi/2, m);
%! [B_comp, D_comp] = ellipticBD(m);
%! assert(abs(B_inc - B_comp) < 1e-13, 'B(pi/2) != B_complete');
%! assert(abs(D_inc - D_comp) < 1e-13, 'D(pi/2) != D_complete');

%!test
%! % Odd symmetry: B(-phi) = -B, D(-phi) = -D, J(-phi) = -J
%! phi = linspace(0.05, 1.2, 30);  m = 0.5;  n = 0.3;
%! [Bp,  Dp,  Jp]  = ellipticBDJ( phi, m, n);
%! [Bm, Dm, Jm] = ellipticBDJ(-phi, m, n);
%! assert(max(abs(Bp + Bm)) < 1e-13, 'B not odd');
%! assert(max(abs(Dp + Dm)) < 1e-13, 'D not odd');
%! assert(max(abs(Jp + Jm)) < 1e-13, 'J not odd');

%!test
%! % Array input: vectorised over phi
%! phi = linspace(0.1, 1.3, 50);  m = 0.7;  n = 0.2;
%! [B, D, J] = ellipticBDJ(phi, m, n);
%! [F, E]    = elliptic12(phi, m .* ones(size(phi)));
%! Pi        = elliptic3(phi, m .* ones(size(phi)), n .* ones(size(phi)));
%! assert(max(abs(B+D-F)) < 1e-12, 'F=B+D array failed');
%! assert(max(abs(B+(1-m)*D-E)) < 1e-12, 'E=B+mc*D array failed');
%! assert(max(abs(B+D+n*J-Pi)) < 1e-10, 'Pi=B+D+nJ array failed');

% ------------------------------------------------------------------
%  Group 2: jacobiEDJ — Jacobi-argument forms
% ------------------------------------------------------------------

%!test
%! % B_u + D_u = u exactly
%! m = 0.5;  n = 0.3;
%! [K, ~] = ellipke(m);
%! u = linspace(0.05, 0.9*K, 20);
%! [Eu, Du] = jacobiEDJ(u, m);
%! Bu = u - Du;
%! assert(max(abs(Bu + Du - u)) < 1e-14, 'B_u+D_u=u failed');

%!test
%! % E_u = u - m*D_u
%! m = 0.5;  n = 0.3;
%! [K, ~] = ellipke(m);
%! u = linspace(0.05, 0.9*K, 20);
%! [Eu, Du] = jacobiEDJ(u, m);
%! assert(max(abs(Eu - (u - m.*Du))) < 1e-14, 'E_u = u-m*D_u failed');

%!test
%! % d(D_u)/du = sn^2(u|m): numerical derivative
%! m = 0.5;  n = 0.3;
%! u0 = 0.6;  h = 1e-7;
%! [~, Du_p] = jacobiEDJ(u0+h, m);
%! [~, Du_m] = jacobiEDJ(u0-h, m);
%! dDu = (Du_p - Du_m) / (2*h);
%! sn = ellipj(u0, m);
%! assert(abs(dDu - sn^2) < 1e-6, ...
%!     sprintf('d(D_u)/du = sn^2 failed: %g vs %g', dDu, sn^2));

%!test
%! % d(E_u)/du = dn^2(u|m): numerical derivative
%! m = 0.5;
%! u0 = 0.6;  h = 1e-7;
%! [Eu_p, ~] = jacobiEDJ(u0+h, m);
%! [Eu_m, ~] = jacobiEDJ(u0-h, m);
%! dEu = (Eu_p - Eu_m) / (2*h);
%! [~, ~, dn] = ellipj(u0, m);
%! assert(abs(dEu - dn^2) < 1e-6, ...
%!     sprintf('d(E_u)/du = dn^2 failed: %g vs %g', dEu, dn^2));

%!test
%! % J_u(u,0|m) = D_u(u|m) [since n=0 makes integrand = sn^2]
%! m = 0.5;
%! [K, ~] = ellipke(m);
%! u = linspace(0.05, 0.9*K, 15);
%! [~, Du, Ju0] = jacobiEDJ(u, m, 0 .* ones(size(u)));
%! assert(max(abs(Ju0 - Du)) < 1e-12, 'J_u(n=0) != D_u');

% ------------------------------------------------------------------
%  Group 3: Carlson RF/RD connection at phi=pi/2
% ------------------------------------------------------------------

%!test
%! % F(pi/2|m) = RF(0, kc^2, 1) = K(m)
%! m = 0.7; kc = sqrt(1-m);
%! [K, ~] = ellipke(m);
%! RF_val = carlsonRF(0, kc^2, 1);
%! assert(abs(RF_val - K) < 1e-13, sprintf('RF(0,kc^2,1) = %g, K = %g', RF_val, K));

%!test
%! % D(pi/2|m) = RD(0, kc^2, 1) / 3 = D(m)
%! m = 0.7; kc = sqrt(1-m);
%! [B_c, D_c] = ellipticBD(m);
%! RD_val = carlsonRD(0, kc^2, 1);
%! assert(abs(RD_val/3 - D_c) < 1e-13, ...
%!     sprintf('RD(0,kc^2,1)/3 = %g, D = %g', RD_val/3, D_c));
