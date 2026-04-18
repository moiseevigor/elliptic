% Tests for Weierstrass elliptic functions.
%
% Test lattice: (e1, e2, e3) = (1, 0, -1)
%   g2 = 4,  g3 = 0,  m = 0.5
%   omega1 = K(0.5)/sqrt(2)  ~  1.31102877714606
%
% Run with:  test testWeierstrass

% ------------------------------------------------------------------
%  Group 1: weierstrassInvariants — round-trip and known values
% ------------------------------------------------------------------

%!test
%! % Invariants for (e1,e2,e3) = (1,0,-1): g2=4, g3=0, Delta=64
%! [g2, g3, D] = weierstrassInvariants(1, 0, -1);
%! assert(abs(g2 - 4)  < 1e-14, 'g2 wrong');
%! assert(abs(g3 - 0)  < 1e-14, 'g3 wrong');
%! assert(abs(D  - 64) < 1e-12, 'Delta wrong');

%!test
%! % Invariants for (0.5, 0, -0.5): g2=1, g3=0
%! [g2, g3, ~] = weierstrassInvariants(0.5, 0, -0.5);
%! assert(abs(g2 - 1) < 1e-14, 'g2 wrong for (0.5,0,-0.5)');
%! assert(abs(g3 - 0) < 1e-14, 'g3 wrong for (0.5,0,-0.5)');

%!test
%! % Array input: invariants computed element-wise
%! e1 = [1, 2]; e2 = [0, 0]; e3 = [-1, -2];
%! [g2, g3, ~] = weierstrassInvariants(e1, e2, e3);
%! assert(isequal(size(g2), [1 2]), 'g2 size wrong');
%! assert(max(abs(g2 - [4, 16])) < 1e-12, 'g2 array wrong');
%! assert(max(abs(g3 - [0,  0])) < 1e-12, 'g3 array wrong');

% ------------------------------------------------------------------
%  Group 2: weierstrassP — special values at half-period (A&S 18.1.4)
% ------------------------------------------------------------------

%!test
%! % P(omega1) = e1 exactly (A&S 18.1.4)
%! e1=1; e2=0; e3=-1;
%! m = (e2-e3)/(e1-e3);
%! [KK,~] = ellipke(m);
%! omega1 = KK / sqrt(e1-e3);
%! P_at_omega1 = weierstrassP(omega1, e1, e2, e3);
%! assert(abs(P_at_omega1 - e1) < 1e-10, ...
%!     sprintf('P(omega1)=%g, expected e1=%g', P_at_omega1, e1));

%!test
%! % P is symmetric: P(z) = P(-z) (even function, A&S 18.2.3)
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.2, 50);
%! Pp = weierstrassP( z, e1, e2, e3);
%! Pm = weierstrassP(-z, e1, e2, e3);
%! assert(max(abs(Pp - Pm)) < 1e-13, 'P is not even');

%!test
%! % P(z) -> Inf as z -> 0 (pole of order 2)
%! e1=1; e2=0; e3=-1;
%! P_small = weierstrassP(1e-10, e1, e2, e3);
%! assert(isinf(P_small) || P_small > 1e15, 'P should blow up near z=0');

%!test
%! % Scalar broadcasting: scalar e values, vector z
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.1, 1.0, 20);
%! P = weierstrassP(z, e1, e2, e3);
%! assert(isequal(size(P), size(z)), 'P size mismatch with scalar e');

%!test
%! % Tensor input: meshgrid z and lattice parameters
%! [zv, e1v] = meshgrid(linspace(0.1, 1.0, 8), [0.5, 1.0, 1.5]);
%! e2v = zeros(size(e1v));
%! e3v = -e1v;
%! P = weierstrassP(zv, e1v, e2v, e3v);
%! assert(isequal(size(P), size(zv)), 'P size mismatch for tensor input');
%! assert(all(isfinite(P(:))), 'P has non-finite values in interior');

% ------------------------------------------------------------------
%  Group 3: weierstrassPPrime — derivative identities
% ------------------------------------------------------------------

%!test
%! % P'(omega1) = 0 exactly (simple zero of ODE RHS at half-period)
%! e1=1; e2=0; e3=-1;
%! m = (e2-e3)/(e1-e3);
%! [KK,~] = ellipke(m);
%! omega1 = KK / sqrt(e1-e3);
%! dP = weierstrassPPrime(omega1, e1, e2, e3);
%! assert(abs(dP) < 1e-10, sprintf('P prime(omega1) = %g, expected 0', dP));

%!test
%! % ODE identity: (P')^2 = 4*P^3 - g2*P - g3  (A&S 18.1.5)
%! e1=1; e2=0; e3=-1;
%! [g2, g3, ~] = weierstrassInvariants(e1, e2, e3);
%! z = linspace(0.05, 1.25, 80);
%! P  = weierstrassP(z, e1, e2, e3);
%! dP = weierstrassPPrime(z, e1, e2, e3);
%! LHS = dP.^2;
%! RHS = 4.*P.^3 - g2.*P - g3;
%! % Exclude points too close to the pole (|P|>1e6)
%! ok = isfinite(P) & abs(P) < 1e6;
%! assert(max(abs(LHS(ok) - RHS(ok))) < 1e-7, ...
%!     sprintf('ODE identity max error: %g', max(abs(LHS(ok)-RHS(ok)))));

%!test
%! % P' is odd: P'(-z) = -P'(z)
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.2, 50);
%! assert(max(abs(weierstrassPPrime( z,e1,e2,e3) + ...
%!             weierstrassPPrime(-z,e1,e2,e3))) < 1e-13, 'P'' not odd');

%!test
%! % Numerical derivative matches weierstrassPPrime
%! e1=1; e2=0; e3=-1;
%! z0 = 0.7;  h = 1e-7;
%! fd = (weierstrassP(z0+h, e1,e2,e3) - weierstrassP(z0-h, e1,e2,e3)) / (2*h);
%! dP = weierstrassPPrime(z0, e1, e2, e3);
%! assert(abs(fd - dP) < 1e-6, ...
%!     sprintf('FD derivative mismatch: fd=%g, PPrime=%g', fd, dP));

% ------------------------------------------------------------------
%  Group 4: weierstrassZeta — derivative and odd-function identities
% ------------------------------------------------------------------

%!test
%! % zeta'(z) = -P(z): numerical derivative check
%! e1=1; e2=0; e3=-1;
%! z0 = 0.6;  h = 1e-6;
%! Zp = (weierstrassZeta(z0+h, e1,e2,e3) - weierstrassZeta(z0-h, e1,e2,e3)) / (2*h);
%! P0 = weierstrassP(z0, e1, e2, e3);
%! assert(abs(Zp + P0) < 1e-5, ...
%!     sprintf('zeta derivative identity failed: Zp=%g, -P=%g', Zp, -P0));

%!test
%! % zeta is odd: zeta(-z) = -zeta(z)
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.1, 30);
%! Zp = weierstrassZeta( z, e1, e2, e3);
%! Zm = weierstrassZeta(-z, e1, e2, e3);
%! assert(max(abs(Zp + Zm)) < 1e-8, 'zeta is not odd');

%!test
%! % zeta(z) -> Inf as z -> 0
%! e1=1; e2=0; e3=-1;
%! Z_small = weierstrassZeta(1e-10, e1, e2, e3);
%! assert(isinf(Z_small) || Z_small > 1e8, 'zeta should blow up near z=0');

%!test
%! % Quasi-periodicity: zeta(z + 2*omega1) = zeta(z) + 2*eta1
%! e1=1; e2=0; e3=-1;
%! m = (e2-e3)/(e1-e3);
%! [KK,~] = ellipke(m);
%! omega1 = KK / sqrt(e1-e3);
%! eta1   = weierstrassZeta(omega1, e1, e2, e3);
%! z0     = 0.5;
%! Z1 = weierstrassZeta(z0,              e1, e2, e3);
%! Z2 = weierstrassZeta(z0 + 2*omega1,   e1, e2, e3);
%! assert(abs(Z2 - Z1 - 2*eta1) < 1e-7, ...
%!     sprintf('quasi-periodicity failed: diff=%g', abs(Z2-Z1-2*eta1)));

% ------------------------------------------------------------------
%  Group 5: weierstrassSigma — logarithmic derivative identity
% ------------------------------------------------------------------

%!test
%! % sigma(0) = 0 exactly
%! e1=1; e2=0; e3=-1;
%! S0 = weierstrassSigma(0, e1, e2, e3);
%! assert(S0 == 0 || abs(S0) < 1e-14, 'sigma(0) should be 0');

%!test
%! % sigma'(0) = 1: numerical derivative check
%! e1=1; e2=0; e3=-1;
%! h  = 1e-7;
%! Sp = (weierstrassSigma(h, e1,e2,e3) - weierstrassSigma(-h, e1,e2,e3)) / (2*h);
%! assert(abs(Sp - 1) < 1e-6, sprintf('sigma prime(0) = %g, expected 1', Sp));

%!test
%! % sigma'/sigma = zeta: logarithmic derivative identity
%! e1=1; e2=0; e3=-1;
%! z0 = 0.6;  h = 1e-6;
%! S0 = weierstrassSigma(z0, e1, e2, e3);
%! Sp = (weierstrassSigma(z0+h, e1,e2,e3) - weierstrassSigma(z0-h, e1,e2,e3)) / (2*h);
%! Z0 = weierstrassZeta(z0, e1, e2, e3);
%! assert(abs(Sp/S0 - Z0) < 1e-5, ...
%!     sprintf('log-deriv identity failed: Sp/S=%g, zeta=%g', Sp/S0, Z0));

%!test
%! % sigma is odd: sigma(-z) = -sigma(z)
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.0, 20);
%! Sp = weierstrassSigma( z, e1, e2, e3);
%! Sm = weierstrassSigma(-z, e1, e2, e3);
%! assert(max(abs(Sp + Sm)) < 1e-8, 'sigma is not odd');

% ------------------------------------------------------------------
%  Group 6: Parallel vs serial consistency
% ------------------------------------------------------------------

%!test
%! % weierstrassP: parallel output matches serial
%! if ~exist('parcellfun'), disp('SKIP: no parallel package'); return; end
%! pkg load parallel;
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.2, 500);
%! elliptic_config('parallel', false);
%! P_s = weierstrassP(z, e1, e2, e3);
%! elliptic_config('parallel', true);
%! P_p = weierstrassP(z, e1, e2, e3);
%! elliptic_config('parallel', false);
%! assert(max(abs(P_s - P_p)) < 1e-13, 'P parallel/serial mismatch');

%!test
%! % weierstrassPPrime: parallel output matches serial
%! if ~exist('parcellfun'), disp('SKIP: no parallel package'); return; end
%! pkg load parallel;
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.2, 500);
%! elliptic_config('parallel', false);
%! dP_s = weierstrassPPrime(z, e1, e2, e3);
%! elliptic_config('parallel', true);
%! dP_p = weierstrassPPrime(z, e1, e2, e3);
%! elliptic_config('parallel', false);
%! assert(max(abs(dP_s - dP_p)) < 1e-13, 'PPrime parallel/serial mismatch');

%!test
%! % weierstrassZeta: parallel output matches serial
%! if ~exist('parcellfun'), disp('SKIP: no parallel package'); return; end
%! pkg load parallel;
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.1, 500);
%! elliptic_config('parallel', false);
%! Z_s = weierstrassZeta(z, e1, e2, e3);
%! elliptic_config('parallel', true);
%! Z_p = weierstrassZeta(z, e1, e2, e3);
%! elliptic_config('parallel', false);
%! assert(max(abs(Z_s - Z_p)) < 1e-10, 'Zeta parallel/serial mismatch');

%!test
%! % weierstrassSigma: parallel output matches serial
%! if ~exist('parcellfun'), disp('SKIP: no parallel package'); return; end
%! pkg load parallel;
%! e1=1; e2=0; e3=-1;
%! z = linspace(0.05, 1.0, 500);
%! elliptic_config('parallel', false);
%! S_s = weierstrassSigma(z, e1, e2, e3);
%! elliptic_config('parallel', true);
%! S_p = weierstrassSigma(z, e1, e2, e3);
%! elliptic_config('parallel', false);
%! assert(max(abs(S_s - S_p)) < 1e-10, 'Sigma parallel/serial mismatch');
