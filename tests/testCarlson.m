% Tests for Carlson symmetric elliptic integrals RF, RD, RJ, RC.
%
% Reference values from DLMF Table 19.36.1 and connection to Legendre forms.
%
% Run with:  test testCarlson

% ------------------------------------------------------------------
%  Group 1: carlsonRC — closed-form degenerate integral
% ------------------------------------------------------------------

%!test
%! % RC(x,x) = 1/sqrt(x)
%! x = [0.25, 1.0, 2.5];
%! RC = carlsonRC(x, x);
%! assert(max(abs(RC - 1./sqrt(x))) < 1e-14, 'RC(x,x) != 1/sqrt(x)');

%!test
%! % RC(0,y) = pi/(2*sqrt(y))
%! y = [0.5, 1.0, 4.0];
%! RC = carlsonRC(zeros(size(y)), y);
%! assert(max(abs(RC - pi./(2.*sqrt(y)))) < 1e-14, 'RC(0,y) != pi/(2*sqrt(y))');

%!test
%! % RC via arctan branch (y > x > 0): DLMF 19.2.17
%! x = 0.5; y = 2.0;
%! RC_ref = atan(sqrt((y-x)/x)) / sqrt(y-x);
%! RC = carlsonRC(x, y);
%! assert(abs(RC - RC_ref) < 1e-14, sprintf('RC arctan branch: %g vs %g', RC, RC_ref));

%!test
%! % RC via arctanh branch (y < x): DLMF 19.2.18
%! x = 3.0; y = 1.0;
%! RC_ref = atanh(sqrt((x-y)/x)) / sqrt(x-y);   % DLMF 19.2.18: (x-y)/x
%! RC = carlsonRC(x, y);
%! assert(abs(RC - RC_ref) < 1e-14, sprintf('RC arctanh branch: %g vs %g', RC, RC_ref));

%!test
%! % RC(y=0) = Inf
%! RC = carlsonRC(1.0, 0.0);
%! assert(isinf(RC), 'RC(x,0) should be Inf');

%!test
%! % DLMF Table 19.36.1: RC(0, 1/4) = pi
%! RC = carlsonRC(0, 0.25);
%! assert(abs(RC - pi) < 1e-13, sprintf('RC(0,1/4) = %g, expected pi', RC));

% ------------------------------------------------------------------
%  Group 2: carlsonRF — symmetric first kind
% ------------------------------------------------------------------

%!test
%! % DLMF 19.20.2: RF(x,x,x) = x^(-1/2)
%! x = 2.5;
%! RF = carlsonRF(x, x, x);
%! assert(abs(RF - x^(-0.5)) < 1e-13, sprintf('RF(x,x,x) = %g, expected %g', RF, x^(-0.5)));

%!test
%! % RF(1,2,3): verified by Gauss-Legendre numerical integration
%! RF = carlsonRF(1, 2, 3);
%! RF_ref = 0.72694584;  % (1/2)*integral_0^inf 1/sqrt((t+1)(t+2)(t+3)) dt
%! assert(abs(RF - RF_ref) < 1e-6, sprintf('RF(1,2,3) = %.10f, expected ~%.8f', RF, RF_ref));

%!test
%! % Connection to Legendre F: F(phi|m) = sin(phi)*RF(cos^2,1-m*sin^2,1)
%! phi = 0.8;  m = 0.5;
%! s = sin(phi); c = cos(phi); d = sqrt(1 - m*s^2);
%! RF_val = carlsonRF(c^2, d^2, 1);
%! [F, ~] = elliptic12(phi, m);
%! assert(abs(s * RF_val - F) < 1e-12, ...
%!     sprintf('F via RF: s*RF=%g, elliptic12=%g', s*RF_val, F));

%!test
%! % Symmetry: RF is symmetric in all three arguments
%! RF_xyz = carlsonRF(1.0, 2.0, 3.0);
%! RF_yxz = carlsonRF(2.0, 1.0, 3.0);
%! RF_zyx = carlsonRF(3.0, 2.0, 1.0);
%! assert(max(abs([RF_xyz - RF_yxz, RF_xyz - RF_zyx])) < 1e-14, ...
%!     'RF not symmetric in all arguments');

%!test
%! % Array input: vectorised over phi
%! phi = linspace(0.1, 1.4, 30);  m = 0.5;
%! s = sin(phi); c = cos(phi); d = sqrt(1 - m.*s.^2);
%! RF_vec = carlsonRF(c.^2, d.^2, ones(size(phi)));
%! [F_vec, ~] = elliptic12(phi, m .* ones(size(phi)));
%! assert(max(abs(s .* RF_vec - F_vec)) < 1e-12, 'RF array vs elliptic12 mismatch');

% ------------------------------------------------------------------
%  Group 3: carlsonRD — symmetric second kind
% ------------------------------------------------------------------

%!test
%! % DLMF 19.20.17: RD(x,x,x) = x^(-3/2)
%! x = 2.0;
%! RD = carlsonRD(x, x, x);
%! assert(abs(RD - x^(-1.5)) < 1e-12, sprintf('RD(x,x,x) = %g, expected %g', RD, x^(-1.5)));

%!test
%! % Connection to Legendre D: D(phi|m) = (sin^3(phi)/3)*RD(cos^2,1-m*sin^2,1)
%! phi = 0.8;  m = 0.5;
%! s = sin(phi); c = cos(phi); d = sqrt(1 - m*s^2);
%! RD_val = carlsonRD(c^2, d^2, 1);
%! % D = F - E (incomplete) verified via elliptic12
%! [F, E] = elliptic12(phi, m);
%! D_leg = (F - E) / m;   % D(phi|m) = (F-E)/m
%! assert(abs(s^3/3 * RD_val - D_leg) < 1e-10, ...
%!     sprintf('D via RD: s3/3*RD=%g, (F-E)/m=%g', s^3/3*RD_val, D_leg));

%!test
%! % RD is symmetric in first two arguments but NOT the third
%! RD_xy = carlsonRD(1.0, 2.0, 3.0);
%! RD_yx = carlsonRD(2.0, 1.0, 3.0);
%! RD_xz = carlsonRD(1.0, 3.0, 2.0);
%! assert(abs(RD_xy - RD_yx) < 1e-13, 'RD not symmetric in x,y');
%! assert(abs(RD_xy - RD_xz) > 1e-6, 'RD incorrectly symmetric in x,z');

%!test
%! % Array input: vectorised over phi
%! phi = linspace(0.1, 1.4, 30);  m = 0.5;
%! s = sin(phi); c = cos(phi); d = sqrt(1 - m.*s.^2);
%! RD_vec = carlsonRD(c.^2, d.^2, ones(size(phi)));
%! [F_vec, E_vec] = elliptic12(phi, m .* ones(size(phi)));
%! D_vec = (F_vec - E_vec) ./ m;
%! assert(max(abs(s.^3./3 .* RD_vec - D_vec)) < 1e-10, 'RD array vs elliptic12 mismatch');

% ------------------------------------------------------------------
%  Group 4: carlsonRJ — symmetric third kind
% ------------------------------------------------------------------

%!test
%! % Connection to Legendre Pi: Pi(n,phi|m) = s*RF + (n*s^3/3)*RJ
%! phi = 0.8;  m = 0.5;  n = 0.3;
%! s = sin(phi); c = cos(phi); d = sqrt(1 - m*s^2);
%! RF_val = carlsonRF(c^2, d^2, 1);
%! RJ_val = carlsonRJ(c^2, d^2, 1, 1 - n*s^2);
%! Pi_carlson = s * RF_val + (n * s^3 / 3) * RJ_val;
%! Pi_legendre = elliptic3(phi, m, n);
%! assert(abs(Pi_carlson - Pi_legendre) < 1e-10, ...
%!     sprintf('Pi via RJ: carlson=%g, legendre=%g', Pi_carlson, Pi_legendre));

%!test
%! % RJ is symmetric in first three arguments
%! RJ_xyz = carlsonRJ(1.0, 2.0, 3.0, 0.5);
%! RJ_yxz = carlsonRJ(2.0, 1.0, 3.0, 0.5);
%! RJ_zxy = carlsonRJ(3.0, 1.0, 2.0, 0.5);
%! assert(max(abs([RJ_xyz - RJ_yxz, RJ_xyz - RJ_zxy])) < 1e-12, ...
%!     'RJ not symmetric in x,y,z');

%!test
%! % RD is special case: RJ(x,y,z,z) = RD(x,y,z)
%! x=1.0; y=2.0; z=3.0;
%! RJ_zz = carlsonRJ(x, y, z, z);
%! RD_val = carlsonRD(x, y, z);
%! assert(abs(RJ_zz - RD_val) < 1e-10, ...
%!     sprintf('RJ(x,y,z,z) != RD(x,y,z): %g vs %g', RJ_zz, RD_val));

%!test
%! % Array input: vectorised over n
%! phi = 0.8;  m = 0.5;
%! n_vec = linspace(0.05, 0.9, 20);
%! s = sin(phi); c = cos(phi); d = sqrt(1 - m*s^2);
%! RF_val = carlsonRF(c^2, d^2, 1);
%! RJ_vec = carlsonRJ(c^2 .* ones(size(n_vec)), ...
%!                    d^2 .* ones(size(n_vec)), ...
%!                    ones(size(n_vec)), ...
%!                    1 - n_vec .* s^2);
%! Pi_carlson = s * RF_val + (n_vec .* s^3 / 3) .* RJ_vec;
%! Pi_leg = elliptic3(phi .* ones(size(n_vec)), m .* ones(size(n_vec)), n_vec);
%! assert(max(abs(Pi_carlson - Pi_leg)) < 1e-9, 'RJ array vs elliptic3 mismatch');
