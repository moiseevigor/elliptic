%TESTEDGECASES  Edge-case tests anchored on values external to this library.
%
% Every expected value below comes from one of three sources that are
% independent of the code under test:
%
%   (1) Closed forms from the reference literature — Abramowitz & Stegun
%       (A&S) and the NIST Digital Library of Mathematical Functions
%       (DLMF, https://dlmf.nist.gov) — evaluated here from GAMMA and PI.
%   (2) Exact functional identities from the same references (quasi-period,
%       Legendre's relation, Jacobi's theta identity, ...).
%   (3) Independent numerical quadrature (QUADGK), which shares no code with
%       any elliptic-integral routine in this library.
%
% The cases chosen are the ones where this library has historically gone
% wrong: phase beyond pi/2, exact multiples of pi/2 (and their neighbouring
% floating-point values), the line Re(u) = pi/2 in the complex plane, and
% m -> 0 / m -> 1.

% ---------------------------------------------------------------------
% A1. K(1/2) and E(1/2) in closed form.
%     K(m=1/2) = Gamma(1/4)^2 / (4*sqrt(pi))                  [A&S 17.3.x]
%     Legendre's relation at k = k' = 1/sqrt(2) gives
%     E(1/2) = K(1/2)/2 + pi/(4*K(1/2))            [A&S 17.3.13, DLMF 19.7.1]
% ---------------------------------------------------------------------
%!test
%! clear
%! K_closed = gamma(0.25)^2 / (4*sqrt(pi));
%! E_closed = K_closed/2 + pi/(4*K_closed);
%! % literature decimals (lemniscate constant omega = 2.6220575542921205,
%! % K(1/2) = omega/sqrt(2)):
%! assert(abs(K_closed - 1.8540746773013723) < 1e-15, 'closed form drifted');
%! assert(abs(E_closed - 1.3506438810476755) < 1e-15, 'closed form drifted');
%! [F,E] = elliptic12(pi/2, 0.5);
%! assert(abs(F - K_closed) < 1e-14, 'F(pi/2|1/2) = %.16g, expected %.16g', F, K_closed);
%! assert(abs(E - E_closed) < 1e-14, 'E(pi/2|1/2) = %.16g, expected %.16g', E, E_closed);
%! [Kk,Ee] = ellipke(0.5);
%! assert(abs(Kk - K_closed) < 1e-14 && abs(Ee - E_closed) < 1e-14, 'ellipke disagrees with closed form');

% ---------------------------------------------------------------------
% A2. Legendre's relation  E*K' + E'*K - K*K' = pi/2   [A&S 17.3.13]
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.05 0.2 0.5 0.8 0.95]
%!     [K,E]   = elliptic12(pi/2, m);
%!     [Kp,Ep] = elliptic12(pi/2, 1-m);
%!     assert(abs(E*Kp + Ep*K - K*Kp - pi/2) < 1e-13, ...
%!         'Legendre relation violated at m=%g by %g', m, E*Kp+Ep*K-K*Kp-pi/2);
%! end

% ---------------------------------------------------------------------
% A3. Degenerate parameters.  m = 0: F = E = phi.
%     m = 1: F(phi|1) = atanh(sin phi), E(phi|1) = sin phi    [A&S 17.4]
% ---------------------------------------------------------------------
%!test
%! clear
%! phi = [-2.5 -0.4 0 0.4 1.2 1.5];
%! [F,E] = elliptic12(phi, zeros(size(phi)));
%! assert(norm(F - phi) < 1e-15 && norm(E - phi) < 1e-15, 'm=0 must give F=E=phi');
%! phi = [-1.2 -0.3 0 0.3 1.2];
%! [F,E] = elliptic12(phi, ones(size(phi)));
%! assert(norm(F - atanh(sin(phi))) < 1e-13, 'F(phi|1) must equal atanh(sin phi)');
%! assert(norm(E - sin(phi))        < 1e-13, 'E(phi|1) must equal sin(phi)');

% ---------------------------------------------------------------------
% B1. Quasi-period.  F(phi + k*pi | m) = F(phi|m) + 2k*K(m),
%     E(phi + k*pi | m) = E(phi|m) + 2k*E(m).                 [A&S 17.4.3]
%     This is the identity that issue #35 violated (the code added k*K).
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.1 0.4 0.75 0.99]
%!     [K,Ec] = elliptic12(pi/2, m);
%!     for phi = [-1.3 -0.2 0.4 1.1 pi/2]
%!         [F0,E0] = elliptic12(phi, m);
%!         for k = -3:3
%!             [Fk,Ek] = elliptic12(phi + k*pi, m);
%!             assert(abs(Fk - (F0 + 2*k*K))  < 1e-12, ...
%!                 'F quasi-period broken: m=%g phi=%g k=%d  got %.15g want %.15g', ...
%!                 m, phi, k, Fk, F0 + 2*k*K);
%!             assert(abs(Ek - (E0 + 2*k*Ec)) < 1e-12, ...
%!                 'E quasi-period broken: m=%g phi=%g k=%d  got %.15g want %.15g', ...
%!                 m, phi, k, Ek, E0 + 2*k*Ec);
%!         end
%!     end
%! end

% ---------------------------------------------------------------------
% B2. Exact multiples of pi/2 and their floating-point neighbours.  These
%     are the knife edges where the period reduction and the Landen branch
%     term pi*ceil(phi/pi-0.5) can jump.  Reference: QUADGK.
% ---------------------------------------------------------------------
%!test
%! clear
%! m = 0.5;
%! pts = [];
%! for j = -5:5
%!     b = j*pi/2;
%!     e = eps(max(abs(b),1));
%!     pts = [pts, b, b+e, b-e, b+3*e, b-3*e];
%! end
%! for u = pts
%!     [F,E] = elliptic12(u, m);
%!     Fq = quadgk(@(t) 1./sqrt(1-m*sin(t).^2), 0, u, 'AbsTol', 1e-13, 'RelTol', 1e-13);
%!     Eq = quadgk(@(t)  sqrt(1-m*sin(t).^2),   0, u, 'AbsTol', 1e-13, 'RelTol', 1e-13);
%!     assert(abs(F-Fq) < 1e-11, 'F(%.17g) = %.15g, quadrature %.15g', u, F, Fq);
%!     assert(abs(E-Eq) < 1e-11, 'E(%.17g) = %.15g, quadrature %.15g', u, E, Eq);
%! end

% ---------------------------------------------------------------------
% B3. Odd symmetry F(-phi) = -F(phi), E(-phi) = -E(phi), including past pi/2.
% ---------------------------------------------------------------------
%!test
%! clear
%! phi = [0.3 1.2 pi/2 2.4 pi 5.0 7.7];
%! for m = [0.2 0.6 0.9]
%!     [Fp,Ep] = elliptic12( phi, m*ones(size(phi)));
%!     [Fm,Em] = elliptic12(-phi, m*ones(size(phi)));
%!     assert(norm(Fp+Fm) < 1e-12 && norm(Ep+Em) < 1e-12, 'parity broken at m=%g', m);
%! end

% ---------------------------------------------------------------------
% C1. Complex phase on the line Re(u) = pi/2 (issue #35).  A&S 17.4.11-13
%     decompose F(phi+i*psi|m) = F(lambda|m) + i*F(mu|1-m); at phi = pi/2
%     and |psi| below the branch point at pi/2 + i*acosh(1/sqrt(m)) this
%     collapses to the closed form
%         F(pi/2 + i*psi | m) = K(m) + i*F(mu | 1-m),
%         sin(mu) = tanh(psi)/sqrt(1-m).
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.1 0.4 0.75 0.9]
%!     psi_c = acosh(1/sqrt(m));                 % branch point ordinate
%!     for frac = [0.005 0.05 0.3 0.9]
%!         psi  = frac*psi_c;
%!         mu   = asin(tanh(psi)/sqrt(1-m));
%!         want = ellipke(m) + 1i*elliptic12(mu, 1-m);
%!         got  = elliptic12i(pi/2 + 1i*psi, m);
%!         assert(abs(got-want) < 1e-12, ...
%!             'F(pi/2+%gi | %g) = %.15g%+.15gi, expected %.15g%+.15gi', ...
%!             psi, m, real(got), imag(got), real(want), imag(want));
%!         % odd in u:  F(-u) = -F(u)
%!         assert(abs(elliptic12i(-(pi/2 + 1i*psi), m) + got) < 1e-12, 'F not odd');
%!     end
%! end

% ---------------------------------------------------------------------
% C2. Continuity of the complex integrals across Re(u) = pi/2 + k*pi while
%     staying below the branch point (|Im u| < acosh(1/sqrt(m))).  A jump
%     of order K(m) here is the visible symptom of issue #35.
% ---------------------------------------------------------------------
%!test
%! clear
%! m = 0.4;
%! psi = 0.5;                                     % acosh(1/sqrt(0.4)) = 1.0317
%! phi = linspace(-4.5, 4.5, 451);
%! [Fi,Ei] = elliptic12i(phi + 1i*psi, m*ones(size(phi)));
%! h = phi(2)-phi(1);
%! assert(max(abs(diff(Fi))) < 5*h, 'F jumps by %g across a half-period', max(abs(diff(Fi))));
%! assert(max(abs(diff(Ei))) < 5*h, 'E jumps by %g across a half-period', max(abs(diff(Ei))));

% ---------------------------------------------------------------------
% C3. On the real axis elliptic12i must reproduce elliptic12 exactly,
%     including past pi/2, and sin(u) must be recovered by sn(F(u|m)|m)
%     (definition of the amplitude, DLMF 22.16).
% ---------------------------------------------------------------------
%!test
%! clear
%! m = 0.4;
%! phi = [-3.4 -1.2 0.3 pi/2 2.0 3.4];
%! [Fr,Er] = elliptic12 (phi, m*ones(size(phi)));
%! [Fc,Ec] = elliptic12i(phi, m*ones(size(phi)));
%! assert(norm(Fc-Fr) < 1e-9 && norm(Ec-Er) < 1e-9, 'elliptic12i /= elliptic12 on the real axis');
%! u = [0.4+0.3i 1.2-0.8i 2.9+0.2i -2.0+0.6i];
%! F = elliptic12i(u, m*ones(size(u)));
%! sn = ellipji(F, m*ones(size(u)));
%! assert(norm(sn - sin(u)) < 1e-9, 'sn(F(u|m)|m) must equal sin(u)');

% ---------------------------------------------------------------------
% D. Jacobi amplitude.  am(0)=0, am(K)=pi/2, am(2K)=pi and
%    am(u + 2K) = am(u) + pi                                  [DLMF 22.16]
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.2 0.6 0.95]
%!     K = ellipke(m);
%!     [~,~,~,a0] = ellipj(0,   m);  assert(abs(a0)        < 1e-14, 'am(0) /= 0');
%!     [~,~,~,a1] = ellipj(K,   m);  assert(abs(a1 - pi/2) < 1e-12, 'am(K) /= pi/2');
%!     [~,~,~,a2] = ellipj(2*K, m);  assert(abs(a2 - pi)   < 1e-12, 'am(2K) /= pi');
%!     for u = [-3.3 -0.7 0.5 2.2 6.0]
%!         [sn,cn,~,am ] = ellipj(u,       m);
%!         [~, ~, ~,am2] = ellipj(u + 2*K, m);
%!         assert(abs(am2 - (am + pi)) < 1e-12, 'am(u+2K) /= am(u)+pi at u=%g m=%g', u, m);
%!         assert(abs(sin(am) - sn) < 1e-13 && abs(cos(am) - cn) < 1e-13, ...
%!             'sin/cos(am) must equal sn/cn at u=%g m=%g', u, m);
%!     end
%! end

% ---------------------------------------------------------------------
% E1. Associate integrals B, D, J.  DLMF 19.25:
%       F   = B + D,  E = B + (1-m)*D,  Pi = B + D + n*J
%     checked past pi/2, and against independent quadrature.
% ---------------------------------------------------------------------
%!test
%! clear
%! n = 0.3;
%! for m = [0.3 0.8]
%!     for phi = [0.7 pi/2 2.0 pi 4.5 2*pi 7.0]
%!         [B,D,J] = ellipticBDJ(phi, m, n);
%!         [F,E]   = elliptic12(phi, m);
%!         assert(abs(B + D - F) < 1e-12, 'B+D /= F at phi=%g m=%g', phi, m);
%!         assert(abs(B + (1-m)*D - E) < 1e-12, 'B+(1-m)D /= E at phi=%g m=%g', phi, m);
%!         Bq = quadgk(@(t) cos(t).^2./sqrt(1-m*sin(t).^2), 0, phi, 'AbsTol',1e-13);
%!         Dq = quadgk(@(t) sin(t).^2./sqrt(1-m*sin(t).^2), 0, phi, 'AbsTol',1e-13);
%!         Jq = quadgk(@(t) sin(t).^2./((1-n*sin(t).^2).*sqrt(1-m*sin(t).^2)), 0, phi, 'AbsTol',1e-13);
%!         assert(abs(B-Bq) < 1e-11, 'B(%g|%g) = %.15g, quadrature %.15g', phi, m, B, Bq);
%!         assert(abs(D-Dq) < 1e-11, 'D(%g|%g) = %.15g, quadrature %.15g', phi, m, D, Dq);
%!         assert(abs(J-Jq) < 1e-11, 'J(%g,%g|%g) = %.15g, quadrature %.15g', phi, n, m, J, Jq);
%!     end
%! end

% ---------------------------------------------------------------------
% E2. Quasi-period of the associate integrals: all three integrands are
%     pi-periodic, so B(phi+k*pi) = B(phi) + 2k*B(m), likewise D and J.
% ---------------------------------------------------------------------
%!test
%! clear
%! m = 0.6; n = 0.25;
%! [Bc,Dc,Jc] = ellipticBDJ(pi/2, m, n);
%! for phi = [-1.1 0.0 0.9 pi/2]
%!     [B0,D0,J0] = ellipticBDJ(phi, m, n);
%!     for k = -2:2
%!         [Bk,Dk,Jk] = ellipticBDJ(phi + k*pi, m, n);
%!         assert(abs(Bk - (B0 + 2*k*Bc)) < 1e-12, 'B quasi-period broken (phi=%g k=%d)', phi, k);
%!         assert(abs(Dk - (D0 + 2*k*Dc)) < 1e-12, 'D quasi-period broken (phi=%g k=%d)', phi, k);
%!         assert(abs(Jk - (J0 + 2*k*Jc)) < 1e-12, 'J quasi-period broken (phi=%g k=%d)', phi, k);
%!     end
%! end

% ---------------------------------------------------------------------
% E3. Jacobi-argument associates.  E_u(u|m) = E(am(u|m) | m) and
%     B_u + D_u = u                                           [DLMF 22.16]
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.3 0.8]
%!     K = ellipke(m);
%!     for u = [0.5 K 1.5*K 2*K 3*K 9.0]
%!         [Eu,Du] = jacobiEDJ(u, m);
%!         [~,~,~,am] = ellipj(u, m);
%!         [~,Eam] = elliptic12(am, m);
%!         assert(abs(Eu - Eam) < 1e-11, 'E_u(%g|%g) /= E(am|m): %.15g vs %.15g', u, m, Eu, Eam);
%!         Dq = quadgk(@(v) arrayfun(@(x) ellipj(x,m)^2, v), 0, u, 'AbsTol', 1e-12);
%!         assert(abs(Du - Dq) < 1e-10, 'D_u(%g|%g) = %.15g, quadrature %.15g', u, m, Du, Dq);
%!     end
%! end

% ---------------------------------------------------------------------
% F. Inverse of the incomplete second kind.  Definitional round trip over
%    several periods and close to m = 1, plus odd symmetry.
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.1 0.5 0.9 0.99]
%!     [~,E1] = ellipke(m);
%!     z  = E1*[-3 -2 -1 -0.5 0 0.01 0.5 0.999 1 1.5 2 3];
%!     iv = inverselliptic2(z, m*ones(size(z)));
%!     [~,Eb] = elliptic12(iv, m*ones(size(z)));
%!     assert(norm(Eb - z, inf) < 1e-11, ...
%!         'inverselliptic2 round trip failed at m=%g by %g', m, norm(Eb-z,inf));
%!     ivm = inverselliptic2(-z, m*ones(size(z)));
%!     assert(norm(ivm + iv, inf) < 1e-11, 'inverselliptic2 must be odd at m=%g', m);
%! end

% ---------------------------------------------------------------------
% G. Complete third kind.  DLMF 19.6:  Pi(0,m) = K(m),
%    Pi(n=m, m) = E(m)/(1-m).
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.1 0.4 0.8]
%!     [K,E] = ellipke(m);
%!     assert(abs(elliptic3(pi/2, m, 0) - K) < 1e-12, 'Pi(0,m) /= K(m) at m=%g', m);
%!     assert(abs(elliptic3(pi/2, m, m) - E/(1-m)) < 1e-11, ...
%!         'Pi(m,m) /= E(m)/(1-m) at m=%g: %.15g vs %.15g', m, elliptic3(pi/2,m,m), E/(1-m));
%! end

% ---------------------------------------------------------------------
% H. Theta functions.  Null values (A&S 16.38):
%      theta2(0)^2 = 2*sqrt(m)*K/pi,  theta3(0)^2 = 2*K/pi,
%      theta4(0)^2 = 2*sqrt(1-m)*K/pi
%    and Jacobi's identity theta3^4 = theta2^4 + theta4^4 (DLMF 20.7.1).
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.05 0.3 0.5 0.9]
%!     K = ellipke(m);
%!     t2 = theta(2,0,m); t3 = theta(3,0,m); t4 = theta(4,0,m);
%!     assert(abs(t2^2 - 2*sqrt(m)  *K/pi) < 1e-13, 'theta2 null value wrong at m=%g', m);
%!     assert(abs(t3^2 - 2         *K/pi) < 1e-13, 'theta3 null value wrong at m=%g', m);
%!     assert(abs(t4^2 - 2*sqrt(1-m)*K/pi) < 1e-13, 'theta4 null value wrong at m=%g', m);
%!     assert(abs(t3^4 - (t2^4 + t4^4)) < 1e-13, 'Jacobi theta identity wrong at m=%g', m);
%!     assert(abs(theta(1,0,m)) < 1e-14, 'theta1(0) must vanish');
%! end

% ---------------------------------------------------------------------
% H2. Jacobi's Theta/Eta versus sn (A&S 16.31):  sn(u|m) = H(u)/(m^(1/4)*Theta(u)),
%     checked past the half period where the old implementation perturbed
%     its own inputs.
% ---------------------------------------------------------------------
%!test
%! clear
%! m = 0.5; K = ellipke(m);
%! for u = [0.3 K 1.5*K 2*K 3*K 5.0 -2.2]
%!     [Th,H] = jacobiThetaEta(u, m);
%!     sn = ellipj(u, m);
%!     assert(abs(H/(m^0.25*Th) - sn) < 1e-12, ...
%!         'H/(m^(1/4) Theta) /= sn at u/K=%g', u/K);
%! end

% ---------------------------------------------------------------------
% I. Carlson symmetric forms.  DLMF 19.20 special values:
%      R_F(x,x,x) = x^(-1/2),  R_D(x,x,x) = x^(-3/2),
%      R_J(x,x,x,x) = x^(-3/2), R_C(x,x) = x^(-1/2),
%      R_F(0,1,1) = R_C(0,1) = pi/2,
%    plus the Legendre connection R_F(0,1-m,1) = K(m).
% ---------------------------------------------------------------------
%!test
%! clear
%! for x = [0.25 1 3.7]
%!     assert(abs(carlsonRF(x,x,x)   - x^-0.5) < 1e-14, 'R_F(x,x,x) wrong at x=%g', x);
%!     assert(abs(carlsonRD(x,x,x)   - x^-1.5) < 1e-14, 'R_D(x,x,x) wrong at x=%g', x);
%!     assert(abs(carlsonRJ(x,x,x,x) - x^-1.5) < 1e-14, 'R_J(x,x,x,x) wrong at x=%g', x);
%!     assert(abs(carlsonRC(x,x)     - x^-0.5) < 1e-14, 'R_C(x,x) wrong at x=%g', x);
%! end
%! assert(abs(carlsonRF(0,1,1) - pi/2) < 1e-14, 'R_F(0,1,1) /= pi/2');
%! assert(abs(carlsonRC(0,1)   - pi/2) < 1e-14, 'R_C(0,1) /= pi/2');
%! for m = [0.1 0.5 0.9]
%!     [K,E] = ellipke(m);
%!     assert(abs(carlsonRF(0,1-m,1) - K) < 1e-13, 'R_F(0,1-m,1) /= K(m) at m=%g', m);
%!     % D(m) = (K-E)/m = R_D(0,1-m,1)/3
%!     assert(abs(carlsonRD(0,1-m,1)/3 - (K-E)/m) < 1e-12, 'R_D(0,1-m,1)/3 /= D(m) at m=%g', m);
%! end

% ---------------------------------------------------------------------
% K. Weierstrass functions against 30-digit reference values.
%    Reference: theta-function closed forms (DLMF 23.6.4, 23.6.8, 23.6.9,
%    23.6.13) evaluated with mpmath 1.4.1 at mp.dps = 50 — a formula
%    family independent of this library's sn-based / series code.  The
%    reference construction itself was validated to < 1e-44 against
%    P'^2 = 4P^3 - g2*P - g3, zeta' = -P, sigma'/sigma = zeta, P(w1) = e1,
%    and the Laurent behaviour at z -> 0 before use.
%    Rows: roots (e1,e2,e3) exactly representable in binary with
%    e1+e2+e3 == 0 exactly (inexact roots break the depressed cubic by
%    O(eps*P^2) and poison the DE identity).
%    z/w1 = 2.6 is the case the old quadrature Sigma got wrong by -277x.
% ---------------------------------------------------------------------
%!test
%! clear
%! % e1 e2 e3, z, P, Pp, Zeta, Sigma  (mpmath, 30 digits)
%! R = [
%!  2.0 -0.5 -1.5 0.364678508599203802930916253782  7.60991654315483213781676575283 -40.716795203467035857794498067   2.73133905217586401231118481172  0.364322893731329436910487536054
%!  2.0 -0.5 -1.5 1.09403552579761140879274876135   2.31139666938461354617370617507   3.65334019657835372774485143298  0.509456906969716181254227726116 0.993298793310168551601976056672
%!  2.0 -0.5 -1.5 2.37041030589482471905095564958   3.55973799836232295109911735041 -11.3205842122331758742412322161   3.5753186226267378852136602119  -7.35675884743112892037462320321
%!  1.0 0.25 -1.25 0.54105576073301731403993176119  3.48955873892538680757168974091 -12.3652547699847269163302027234   1.83475002192072477467486916472  0.540061003157461700899655935519
%!  1.0 0.25 -1.25 1.62316728219905194211979528357  1.13305294818566893955539247187   1.05828457240220452238075426776  0.25696621923989724678288838394  1.40818598628919804944444268506
%!  1.0 0.25 -1.25 3.51686244476461254125955644774  1.67783096402393224759843061699  -3.36668198668738980119457645764  2.26609152750196941530142216878 -8.29563954795637078192140811817
%! ];
%! for i = 1:rows(R)
%!     e1=R(i,1); e2=R(i,2); e3=R(i,3); z=R(i,4);
%!     assert(abs(weierstrassP(z,e1,e2,e3)      - R(i,5)) < 1e-12*max(1,abs(R(i,5))), 'P wrong at row %d', i);
%!     assert(abs(weierstrassPPrime(z,e1,e2,e3) - R(i,6)) < 1e-12*max(1,abs(R(i,6))), 'Pprime wrong at row %d', i);
%!     assert(abs(weierstrassZeta(z,e1,e2,e3)   - R(i,7)) < 1e-12*max(1,abs(R(i,7))), 'Zeta wrong at row %d', i);
%!     assert(abs(weierstrassSigma(z,e1,e2,e3)  - R(i,8)) < 1e-12*max(1,abs(R(i,8))), 'Sigma wrong at row %d', i);
%! end

% ---------------------------------------------------------------------
% K2. Weierstrass structural identities (DLMF 23.2, 23.6):
%     P'^2 = 4P^3 - g2 P - g3;  zeta(z + 2k*w1) = zeta(z) + 2k*eta1 with
%     eta1 = zeta(w1);  sigma(z + 2*w1) = -sigma(z)*exp(2*eta1*(z + w1)).
%     Roots exactly binary-representable, e1+e2+e3 == 0 exactly.
% ---------------------------------------------------------------------
%!test
%! clear
%! e1=1.5; e2=0.25; e3=-1.75;
%! [g2,g3] = weierstrassInvariants(e1,e2,e3);
%! m = (e2-e3)/(e1-e3);  w1 = ellipke(m)/sqrt(e1-e3);
%! eta1 = weierstrassZeta(w1,e1,e2,e3);
%! for zf = [0.23 0.61 0.89]
%!     z = zf*w1;
%!     P  = weierstrassP(z,e1,e2,e3);
%!     Pp = weierstrassPPrime(z,e1,e2,e3);
%!     assert(abs(Pp^2 - (4*P^3 - g2*P - g3)) < 1e-9*max(1,Pp^2), 'DE violated at z/w1=%g', zf);
%!     for k = [-2 -1 1 2]
%!         Zk = weierstrassZeta(z + 2*k*w1, e1,e2,e3);
%!         assert(abs(Zk - (weierstrassZeta(z,e1,e2,e3) + 2*k*eta1)) < 1e-11, ...
%!             'zeta quasi-period broken at z/w1=%g k=%d', zf, k);
%!     end
%!     Sq = weierstrassSigma(z + 2*w1, e1,e2,e3);
%!     Sw = -weierstrassSigma(z,e1,e2,e3) * exp(2*eta1*(z + w1));
%!     assert(abs(Sq - Sw) < 1e-10*max(1,abs(Sw)), 'sigma quasi-period broken at z/w1=%g', zf);
%! end

% ---------------------------------------------------------------------
% J. Ellipse arclength.  Degenerate circle and the full ellipse, against
%    the plain arclength integral (QUADGK) — the full ellipse spans 2*pi
%    of phase and so exercises the quasi-period.
% ---------------------------------------------------------------------
%!test
%! clear
%! assert(abs(arclength_ellipse(2,2) - 2*pi*2) < 1e-12, 'circle perimeter wrong');
%! for ab = [5 10; 10 5; 1 0.5; 1 2]'
%!     a = ab(1); b = ab(2);
%!     q = quadgk(@(t) sqrt(a^2*sin(t).^2 + b^2*cos(t).^2), 0, 2*pi, 'AbsTol', 1e-12);
%!     assert(abs(arclength_ellipse(a,b) - q) < 1e-9, ...
%!         'full ellipse a=%g b=%g: %.15g vs quadrature %.15g', a, b, arclength_ellipse(a,b), q);
%!     for t1 = [pi/10 1.9 pi 4.0]
%!         q1 = quadgk(@(t) sqrt(a^2*sin(t).^2 + b^2*cos(t).^2), 0, t1, 'AbsTol', 1e-12);
%!         assert(abs(arclength_ellipse(a,b,0,t1) - q1) < 1e-9, ...
%!             'arc a=%g b=%g to %g: %.15g vs %.15g', a, b, t1, arclength_ellipse(a,b,0,t1), q1);
%!     end
%! end

% ---------------------------------------------------------------------
% L. Third kind past pi/2 and negative phase — mpmath 1.4.1 ellippi at
%    mp.dps = 50 (values to 25 digits).  The Gauss-Legendre core is only
%    valid on [0, pi/2]; these rows exercise the quasi-period reduction
%    Pi(u+k*pi|m,c) = Pi(u|m,c) + 2k*Pi(pi/2|m,c) and oddness.
% ---------------------------------------------------------------------
%!test
%! clear
%! R = [ % u, c, m, Pi   (mpmath ellippi(c, u, m))
%!  2.0 0.3 0.4 2.906928215899509476566159
%!  3.14159265358979323846 0.3 0.4 4.297590831732913271869402
%!  6.0 0.3 0.4 8.30819635557621288613689
%!  -2.2 0.25 0.6 -3.471643753410082592812564
%!  4.5 0.7 0.2 8.475833763417048864937332
%! ];
%! for i = 1:rows(R)
%!     got = elliptic3(R(i,1), R(i,3), R(i,2));
%!     assert(abs(got - R(i,4)) < 1e-12*max(1,abs(R(i,4))), ...
%!         'Pi(%g|m=%g,c=%g) = %.15g, expected %.15g', R(i,1), R(i,3), R(i,2), got, R(i,4));
%! end
%! % oddness on the extended domain
%! assert(abs(elliptic3(-2.0, 0.4, 0.3) + elliptic3(2.0, 0.4, 0.3)) < 1e-13, 'Pi must be odd');

% ---------------------------------------------------------------------
% M. Complex Jacobi functions — mpmath 1.4.1 ellipfun at mp.dps = 50
%    (values to 25 digits), including multi-period complex arguments.
% ---------------------------------------------------------------------
%!test
%! clear
%! R = [ % m, Re u, Im u, Re sn, Im sn, Re cn, Im cn, Re dn, Im dn
%!  0.3 1.2 -0.8 1.144487529338893062450354 -0.2802393491508895255196492 0.4746455171789828435720278 0.6757262603879129546403857 0.8030935803341932691444772 0.1198106104396107763170419
%!  0.7 2.9 0.2 0.8945661685904884966950071 -0.0623052446149157017922886 -0.4667874916734775923001321 -0.1194037221486780885477261 0.667799853302213119590502 0.05842366478197534982879624
%!  0.7 0.1 1.9 2.711675412858391500357814 -4.920359160836619695635427 -4.998252993291208495782488 -2.669416089337967774258619 -4.209736402055241155323296 -2.218593037476527723830424
%!  0.96 4.5 -2.5 0.9833643701168131893646114 -0.09200718087257114311996756 0.3369732007632264400404546 0.2684978605421883948528415 -0.3680824599803019014518416 -0.2359729940161469769037399
%! ];
%! for i = 1:rows(R)
%!     [sn,cn,dn] = ellipji(R(i,2)+1i*R(i,3), R(i,1));
%!     assert(abs(sn - (R(i,4)+1i*R(i,5))) < 1e-12*max(1,abs(sn)), 'sn wrong at row %d', i);
%!     assert(abs(cn - (R(i,6)+1i*R(i,7))) < 1e-12*max(1,abs(cn)), 'cn wrong at row %d', i);
%!     assert(abs(dn - (R(i,8)+1i*R(i,9))) < 1e-12*max(1,abs(dn)), 'dn wrong at row %d', i);
%! end

% ---------------------------------------------------------------------
% N. Nome — closed form q(1/2) = exp(-pi)  (K(1/2) = K'(1/2)); AGM vs
%    mpmath.agm.  theta_prime spot value against mpmath jtheta.
% ---------------------------------------------------------------------
%!test
%! clear
%! assert(abs(nomeq(0.5) - exp(-pi)) < 1e-16, 'q(1/2) must equal exp(-pi)');
%! assert(abs(inversenomeq(exp(-pi)) - 0.5) < 1e-12, 'inversenomeq(exp(-pi)) must be 1/2');
%! % AGM(1, 0.01) = 0.2621668872022492366947771  (mpmath 1.4.1, dps=50)
%! [a,b,c,n] = agm(1, 0.01, (1-0.01)/2);
%! assert(abs(a(double(n)+1,1) - 0.2621668872022492366947771) < 1e-15, 'AGM(1,0.01) wrong');
%! % theta_1'(0.4 | m=0.6) = jtheta(1, 0.4, q(0.6), 1)  (mpmath, 25 digits)
%! [th, thp] = theta_prime(1, 0.4, 0.6);
%! assert(abs(th  - 0.3776251831225481533690215) < 1e-13, 'theta1(0.4|0.6) wrong');
%! assert(abs(thp - 0.8967180488866961938686533) < 1e-13, 'theta1_prime(0.4|0.6) wrong');

% ---------------------------------------------------------------------
% O. Theta derivative at the zeros of theta itself.  DLMF 20.4.6:
%    theta1'(0) = theta2(0)*theta3(0)*theta4(0).  The old logarithmic-
%    derivative form returned NaN (0*Inf) at z = k*pi for theta1 and
%    z = pi/2 + k*pi for theta2.
% ---------------------------------------------------------------------
%!test
%! clear
%! for m = [0.1 0.5 0.9]
%!     [~, d1] = theta_prime(1, 0, m);
%!     [t2, ~] = theta_prime(2, 0, m);
%!     [t3, ~] = theta_prime(3, 0, m);
%!     [t4, ~] = theta_prime(4, 0, m);
%!     assert(abs(d1 - t2*t3*t4) < 1e-13, 'DLMF 20.4.6 violated at m=%g', m);
%!     [~, dpi] = theta_prime(1, pi, m);
%!     assert(~isnan(dpi) && abs(dpi + d1) < 1e-13, 'theta1''(pi) must equal -theta1''(0)');
%!     [~, d2] = theta_prime(2, pi/2, m);
%!     assert(~isnan(d2), 'theta2''(pi/2) must be finite');
%! end

% ---------------------------------------------------------------------
% P. uniquetol_compat — compatibility utility (the numerical kernels now
%    group exact duplicates only). Contract: C(ic) reconstructs A within tol,
%    C == A(ia) exactly, C strictly increasing, and near-duplicates
%    (within tol) collapse to one group.  A wrong index map here would
%    silently corrupt every m-grouped elliptic value downstream.
% ---------------------------------------------------------------------
%!test
%! clear
%! rand('seed', 42);
%! base = rand(1, 20);
%! A = [base, base + 1e-13, 0.5, 0.5+2e-12, 0.7, fliplr(base)];   % dups + jitter
%! tol = 1e-11;
%! [C, ia, ic] = uniquetol_compat(A, tol);
%! assert(all(abs(C(ic) - A) <= tol * max(1, abs(A))), 'C(ic) must reconstruct A within tol');
%! assert(isequal(C(:), A(ia)(:)), 'C must equal A(ia) exactly');
%! assert(all(diff(C) > 0), 'C must be strictly increasing');
%! assert(all(diff(C) > tol * max(1, abs(C(1:end-1)))), 'groups closer than tol survived');
%! % elliptic12 must preserve every distinct m rather than tolerance-group it:
%! m_dup = [0.4, 0.4+5e-13, 0.4-5e-13, 0.7, 0.7+1e-12];
%! [F1,E1] = elliptic12(1.1*ones(size(m_dup)), m_dup);
%! for k = 1:numel(m_dup)
%!     [F2,E2] = elliptic12(1.1, m_dup(k));
%!     assert(abs(F1(k)-F2) < 2e-13 && abs(E1(k)-E2) < 2e-13, ...
%!         'grouped vs scalar elliptic12 disagree at k=%d', k);
%! end

% ---------------------------------------------------------------------
% Q. Adversarial-review round (external Codex + mpmath 1.4.1, dps=40).
%    Each block below is a counterexample that a prior version failed.
% ---------------------------------------------------------------------
%!test
%! clear
%! % Q1: negative amplitude while the COMPLETE integral has a pole (0*Inf).
%! assert(abs(elliptic3(-1, 0.5, 1)   - (-1.7319915420235269928)) < 1e-13, 'Pi(-1|.5,1) wrong');
%! assert(abs(elliptic3(-1, 1,   0.2) - (-1.3115010674599590753)) < 1e-13, 'Pi(-1|1,.2) wrong');
%! assert(~isnan(elliptic3(-0.4, 1, 1)), 'Pi(-0.4|1,1) must not be NaN');

%!test
%! clear
%! % Q2: complex F/E small-m series region (A&S path lost sqrt(eps/m) digits).
%! assert(abs(elliptic12i(0.2i, 1e-20) - 0.2i) < 1e-15, 'F(0.2i|1e-20)');
%! [F,E] = elliptic12i(pi/2 + 0.2i, 1e-14);
%! assert(abs(F - (1.5707963267949005462 + 0.20000000000000101344i)) < 1e-13, 'F(pi/2+0.2i|1e-14)');
%! assert(abs(E - (1.5707963267948926922 + 0.19999999999999898656i)) < 1e-13, 'E(pi/2+0.2i|1e-14)');
%! [F,E] = elliptic12i(pi/2 + 0.2i, 1e-6);
%! assert(abs(F - (1.5707967194941992113 + 0.20000010134411776594i)) < 5e-12, 'F(pi/2+0.2i|1e-6)');
%! assert(abs(E - (1.5707959340957412894 + 0.19999989865593359446i)) < 5e-12, 'E(pi/2+0.2i|1e-6)');
%! % both sides of the series threshold vs mpmath (dps=30)
%! Fa = elliptic12i(1.1+0.3i, 0.99e-4);
%! assert(abs(Fa - (1.1000153646885162+0.3000120622623928i)) < 1e-12, 'series side of crossover');
%! Fb = elliptic12i(1.1+0.3i, 1.01e-4);
%! assert(abs(Fb - (1.1000156750952820+0.3000123059589823i)) < 5e-12, 'decomposition side of crossover');

%!test
%! clear
%! % Q3: Weierstrass near-origin values are huge but FINITE (DLMF 23.9.2);
%! % only the exact lattice point is a pole.
%! assert(abs(weierstrassP(1e-16,1,0,-1)      - 1e32)  < 1e19,  'P(1e-16) must be ~1e32');
%! assert(abs(weierstrassPPrime(1e-16,1,0,-1) + 2e48)  < 1e36,  'Pp(1e-16) must be ~-2e48');
%! assert(abs(weierstrassZeta(1e-16,1,0,-1)   - 1e16)  < 1e3,   'zeta(1e-16) must be ~1e16');
%! assert(isinf(weierstrassP(0,1,0,-1)), 'P(0) must be Inf');

%!test
%! clear
%! % Q4: inverse nome via DLMF 20.9.1 -- exact at every scale.
%! assert(abs(inversenomeq(1e-30) - 1.6e-29) < 1e-41, 'm(1e-30) must be 1.6e-29');
%! assert(abs(inversenomeq(1e-12) - 1.5999999999872e-11) < 1e-24, 'm(1e-12) = 16q - 128q^2');
%! for mv = [1e-8 0.3 0.85 0.999]
%!     assert(abs(inversenomeq(nomeq(mv)) - mv) < 1e-12*max(mv,1e-3), 'roundtrip fails at m=%g', mv);
%! end

%!test
%! clear
%! % Q5: Carlson scale invariance (DLMF 19.20): RF ~ lambda^(-1/2), RC same,
%! % RD/RJ ~ lambda^(-3/2).  An absolute branch tolerance broke this.
%! x=1; y=2; z=3; p=4;
%! for lam = [1e-20 1e20]
%!     assert(abs(carlsonRF(lam*x,lam*y,lam*z)     - carlsonRF(x,y,z)/sqrt(lam))       < 1e-10*abs(carlsonRF(x,y,z)/sqrt(lam)), 'RF homogeneity at %g', lam);
%!     assert(abs(carlsonRC(lam*x,lam*y)           - carlsonRC(x,y)/sqrt(lam))         < 1e-10*abs(carlsonRC(x,y)/sqrt(lam)), 'RC homogeneity at %g', lam);
%!     assert(abs(carlsonRD(lam*x,lam*y,lam*z)     - carlsonRD(x,y,z)/lam^1.5)         < 1e-10*abs(carlsonRD(x,y,z)/lam^1.5), 'RD homogeneity at %g', lam);
%!     assert(abs(carlsonRJ(lam*x,lam*y,lam*z,lam*p) - carlsonRJ(x,y,z,p)/lam^1.5)     < 1e-10*abs(carlsonRJ(x,y,z,p)/lam^1.5), 'RJ homogeneity at %g', lam);
%! end
%! assert(abs(carlsonRC(1e-20,2e-20) - 7853981633.9744830962) < 1e-4, 'RC(1e-20,2e-20)');

%!test
%! clear
%! % Q6: ellipticBD nondegenerate anchors (mpmath: B=(E-(1-m)K)/m, D=(K-E)/m).
%! R = [0.2 0.8066808960371526438 0.85294270257337535705
%!      0.7 0.88437375336868858245 1.1909893819237805614
%!      0.999 0.99832798626015502386 3.8428045742901420065];
%! for i = 1:rows(R)
%!     [B,D] = ellipticBD(R(i,1));
%!     assert(abs(B - R(i,2)) < 1e-14, 'B(%g)', R(i,1));
%!     assert(abs(D - R(i,3)) < 1e-13, 'D(%g)', R(i,1));
%! end

%!test
%! clear
%! % Q7: reversed arc intervals are signed, circles included.
%! assert(abs(arclength_ellipse(2,3,1,0.1) + arclength_ellipse(2,3,0.1,1)) < 1e-13, 'ellipse arc not odd under reversal');
%! assert(abs(arclength_ellipse(2,2,1,0.1) - (-1.8)) < 1e-13, 'reversed circle arc must be -a*(t1-t0)');

% ---------------------------------------------------------------------
% R. Second adversarial round (fuzz vs mpmath dps=40 over parameter
%    endpoints, extreme scales, poles and period multiples).  Every
%    reference was evaluated at the EXACT DOUBLE the library receives --
%    near singularities the decimal input differs from its double by
%    enough to move the answer at 1e-9.  Each block is a counterexample a
%    prior version failed; scipy reaches machine precision on all of them.
% ---------------------------------------------------------------------
%!test
%! clear
%! m1 = 1 - eps/2;
%! % R1: m -> 1 near phi = pi/2 (Delta^2 formed as (1-m) + m cos^2)
%! assert(abs(elliptic12(pi/2-1e-9, m1) - 19.65993026560449767) < 1e-12*20, 'F(pi/2-1e-9 | 1-eps/2)');
%! % R2: m = 1 exactly: F must be exactly 0 at 0 and odd
%! assert(elliptic12(0, 1) == 0, 'F(0|1) must be exactly 0');
%! assert(abs(elliptic12(1e-16, 1) - 1e-16) < 1e-31, 'F(1e-16|1) must be 1e-16');
%! % R3: third kind at the endpoint poles
%! assert(abs(elliptic3(pi/2-1e-6, 0.3, 1) - 1195228.2584444625825) < 1e-12*1.2e6, 'Pi(pi/2-1e-6 | .3, c=1)');
%! assert(abs(elliptic3(pi/2-1e-6, 1-1e-8, 0.9) - 88.615055050793590585) < 1e-12*90, 'Pi(pi/2-1e-6 | 1-1e-8, .9)');

%!test
%! clear
%! % R4: Carlson -- disparate scales, tiny y, near-equal args, double zeros
%! assert(abs(carlsonRJ(1e-20,2e-20,3e-20,0.5) - 43616756114.805842986) < 1e-12*4.4e10, 'RJ disparate scales');
%! assert(abs(carlsonRJ(2,3,4,1e-10) - 7.179193296087372323) < 1e-13*7.2, 'RJ small p');
%! assert(abs(carlsonRC(3,1e-10) - 7.3643213780616827229) < 1e-13*7.4, 'RC(3,1e-10)');
%! assert(abs(carlsonRC(1.0000000000001,1) - 0.99999999999998334665) < 1e-14, 'RC(1+1e-13,1)');
%! assert(isinf(carlsonRF(0,0,1)) && isinf(carlsonRD(0,0,1)) && isinf(carlsonRJ(0,0,1,2)), 'two zero args must be Inf');

%!test
%! clear
%! % R5: Jacobi functions at m -> 1 (atan2 Landen step)
%! assert(abs(ellipj(9.375277798108883, 1-eps/2) - 0) >= 0);   % smoke: callable
%! [~,cn] = ellipj(9.375277798108883, 1-eps/2);
%! assert(abs(cn - 0.00016958935096417269446) < 1e-13*1.7e-4, 'cn(9.375 | 1-eps/2)');
%! [~,cn] = ellipj(7, 1-1e-12);
%! assert(abs(cn - 0.0018237622775256289351) < 1e-13*1.8e-3, 'cn(7 | 1-1e-12)');

%!test
%! clear
%! % R6: inverse E for tiny negative z (oddness first, relative stop)
%! for m = [0 0.5 1-1e-8]
%!     [~,E1] = ellipke(m);  z = -1e-9*E1;
%!     [~,Eb] = elliptic12(inverselliptic2(z, m), m);
%!     assert(abs(Eb - z) < 1e-13*abs(z), 'inverselliptic2 tiny negative z at m=%g', m);
%! end

%!test
%! clear
%! % R7: complex F -- cancellation-free tan^2(mu): tiny psi and small m
%! assert(abs(imag(elliptic12i(pi/2 + 1e-9i, 0.9)) - 3.1622776601683798848e-9) < 1e-12*3.2e-9, 'Im F(pi/2+1e-9i | .9)');
%! F = elliptic12i(pi/2 + 1e-9i, 1-eps/2);
%! assert(abs(F - (19.754694640120759063 + 0.095049319491958534055i)) < 1e-9*20, 'F(pi/2+1e-9i | 1-eps/2)');
%! F = elliptic12i(0.4 + 0.3i, 1e-4);
%! assert(abs(F - (0.39999936996865927219 + 0.30000195549059365219i)) < 1e-13, 'F(0.4+0.3i | 1e-4)');

%!test
%! clear
%! % R8: nome at tiny m (K' from the exact argument m); Weierstrass with the
%! % period computed from 1-m = (e1-e2)/(e1-e3) on a near-m=1 lattice
%! assert(abs(nomeq(1e-16) - 6.2500000000000001819e-18) < 1e-13*6.25e-18, 'q(1e-16)');   % pi*K'/K ~ 40 eps amplification
%! assert(abs(nomeq(1e-17) - 6.25e-19) < 1e-13*6.25e-19, 'q(1e-17)');
%! e1 = 0.5000001; e2 = 0.5; e3 = -1.0000001;  z = 13.391953465243201;
%! assert(abs(weierstrassP(z,e1,e2,e3)     - 0.51848214450600943279)     < 1e-12, 'P near-1 lattice');
%! assert(abs(weierstrassZeta(z,e1,e2,e3)  - -5.4787546526901492279)  < 1e-12*5.5, 'Zeta near-1 lattice');
%! assert(abs(weierstrassSigma(z,e1,e2,e3) - 1.822626274365935705e-13) < 1e-12*1.8e-13, 'Sigma near-1 lattice');

% ---------------------------------------------------------------------
% S. Third adversarial round: dense random fuzz + API abuse.
%    S1 Cody-Waite pi reduction: F/E/Z/Pi at u = 5.5*pi + 8e-10, m = 1-1.5e-13
%       (mpmath at the exact double inputs).  Before the split, k*pi rounding
%       cost eps*|u| in the reduced phase, amplified 1e5x by dZ/dphi there.
%    S2 NaN in, NaN out (a NaN m used to crash elliptic12's grouping).
%    S3 elliptic12i(-0) is -0 exactly (the cot nudge returned eps).
%    S4 R_J at an argument ratio of 1.9e44 (fixed 100 duplications).
% ---------------------------------------------------------------------
%!test
%! clear
%! u = 17.27875959554386;  m = 0.99999999999985;
%! [F,E,Z] = elliptic12(u, m);
%! assert(abs(F - 177.65640489133312311) < 2e-9,  'F at 5.5pi+8e-10, m->1 (conditioning floor ~1e-9)');
%! assert(abs(E - 11.000000000012911122) < 1e-12, 'E at 5.5pi+8e-10');
%! assert(abs(Z - (-0.00012790056416609388015)) < 1e-10, 'Z at 5.5pi+8e-10: k*pi rounding used to cost 1e-9');
%! assert(abs(elliptic3(u, m, 0.3) - 248.50046674013002377) < 3e-9, 'Pi at 5.5pi+8e-10');

%!test
%! clear
%! assert(isnan(elliptic12(0.3, NaN)) && isnan(elliptic12(NaN, 0.5)), 'elliptic12 must propagate NaN');
%! v = elliptic12([0.3 0.5 0.7], [0.2 NaN 0.4]);
%! assert(isnan(v(2)) && ~any(isnan(v([1 3]))), 'NaN must not leak into neighbours');
%! assert(isnan(ellipj(0.3, NaN)), 'ellipj must propagate NaN');
%! assert(elliptic12i(-0, 0.5) == 0 && elliptic12i(0, 0.5) == 0, 'F(0) must be exactly 0');
%! assert(abs(carlsonRJ(2.798e-18, 5.954e-24, 9.634e-23, 1.134e21) - 9.9678905686736778972e-12) < 1e-12*1e-11, 'RJ at ratio 1.9e44');

%% ---------------------------------------------------------------------
%% T. Theta at a huge argument (mpmath jtheta at the exact double v).
%%    Forming (2n+1)*v as a double product rounded by eps*|k v| (2e-10 here,
%%    9e-9 at v ~ 1e11); the series now uses the angle-addition recurrence,
%%    and theta() no longer round-trips v -> u -> v through jacobiThetaEta.
%% ---------------------------------------------------------------------
%!test
%! clear
%! v = 123456789.123;
%! [t, tp] = theta_prime(1, v, 0.4);
%! assert(abs(t - (0.84585020823346348431)) < 2e-15, 'theta1 at v=1.2e8');
%! assert(abs(tp - (0.015114923736622936955)) < 2e-14, 'theta1'' at v=1.2e8');
%! assert(abs(theta(1, v, 0.4) - (0.84585020823346348431)) < 2e-15, 'theta() at v=1.2e8');
%! [t, tp] = theta_prime(2, v, 0.4);
%! assert(abs(t - (0.014932290326334898348)) < 2e-15, 'theta2 at v=1.2e8');
%! assert(abs(tp - (-0.84241862816186020729)) < 2e-14, 'theta2'' at v=1.2e8');
%! assert(abs(theta(2, v, 0.4) - (0.014932290326334898348)) < 2e-15, 'theta() at v=1.2e8');
%! [t, tp] = theta_prime(3, v, 0.4);
%! assert(abs(t - (0.93627542467710214194)) < 2e-15, 'theta3 at v=1.2e8');
%! assert(abs(tp - (-0.004519191759268069986)) < 2e-14, 'theta3'' at v=1.2e8');
%! assert(abs(theta(3, v, 0.4) - (0.93627542467710214194)) < 2e-15, 'theta() at v=1.2e8');
%! [t, tp] = theta_prime(4, v, 0.4);
%! assert(abs(t - (1.0637286984176921296)) < 2e-15, 'theta4 at v=1.2e8');
%! assert(abs(tp - (0.0045203629452149075654)) < 2e-14, 'theta4'' at v=1.2e8');
%! assert(abs(theta(4, v, 0.4) - (1.0637286984176921296)) < 2e-15, 'theta() at v=1.2e8');

%% ---------------------------------------------------------------------
%% U. Round 6 (cross-port parity sweep at extreme m, 2026-09-02).  Every
%%    anchor is mpmath at the EXACT double inputs.
%%    - m in [eps^2, ~5e-16]: the AGM converges in one step, no Landen step
%%      ran and the scale e stayed 0 -> F = E = Inf.
%%    - theta nome from ellipke(1-m): 1-m rounds, q was 30% off at m ~ 1e-16.
%%    - E/K sum stopped one AGM term early: E off by 1.4e-13 near m -> 1.
%%    - k*pi reduction: k*double(pi) rounds by eps*|u| (2e-10 at u = 1e6);
%%      sub_kpi splits pi into 25-bit parts so k*PI_A, k*PI_B are exact.
%%    - elliptic3 reflection used Pi(double(pi/2)) as the complete integral;
%%      cos(double(pi/2)) = 6e-17 is not 0 and the sliver is 2e-7 at m = 1-eps/2.
%%    - carlsonRJ series: E3 = XYZ + 2 E2 P + 4 P^3 (DLMF 19.36.2), not 3 P^3.
%%    - ellipticBDJ: n > 1 beyond the pole gave complex J silently; n = 1 gave NaN.
%%    - elliptic3 now accepts c < 0 (as the Python port does).
%% ---------------------------------------------------------------------
%!test
%! clear
%! [F, E] = elliptic12([1 1], [3e-16 5e-16]);
%! assert(all(isfinite(F)) && all(isfinite(E)), 'F, E must be finite for m ~ 3e-16 (were Inf)');
%! assert(abs(F(1) - 1.0) < 5e-16 && abs(E(1) - 0.99999999999999996) < 5e-16, 'F(1|3e-16), E(1|3e-16)');
%! assert(abs(F(2) - 1.0000000000000001) < 5e-16 && abs(E(2) - 0.99999999999999993) < 5e-16, 'F(1|5e-16), E(1|5e-16)');
%! assert(abs(theta(1, 34401.9, 1.6e-16) - 0.00011178415088289534) < 1e-13 * 1.1e-4, 'theta1 at m = 1.6e-16 (nome from exact m)');
%! assert(abs(theta_prime(1, 6577.39, 1.5e-16) - (-9.8878892558450512e-5)) < 1e-13 * 1e-4, 'theta_prime at m = 1.5e-16');
%! [~, H] = jacobiThetaEta(6577.39 * 2 * ellipke(1.5e-16) / pi, 1.5e-16);
%! assert(abs(H - (-9.8878892558450512e-5)) < 1e-12 * 1e-4, 'jacobiThetaEta eta at m = 1.5e-16');
%! [F, E] = elliptic12(-1.65181, 0.99999999999999578);
%! assert(abs(E - (-1.0032798131910099)) < 2e-14, 'E near m -> 1 (missing AGM term gave 1.4e-13)');
%! assert(abs(F - (-32.666065762173088)) < 1e-14 * 33, 'F near m -> 1');
%! u = 80101.48788857895;  m = 0.9999533239086507;         % 25497*pi + 0.3
%! [F, E, Z] = elliptic12(u, m);
%! assert(abs(Z - 0.2477141143165845) < 2e-14, 'Jacobi Zeta at u = 8e4 (k*pi split)');
%! assert(abs(E - 51001.284415600044) < 1e-14 * 51001, 'E at u = 8e4');
%! assert(abs(F - 324959.38078465716) < 1e-14 * 324959, 'F at u = 8e4');
%! [~, ~, Z] = elliptic12(1000000.123, 1 - eps/2);
%! assert(abs(Z - (-0.220434859492317)) < 2e-14, 'Jacobi Zeta at u = 1e6, m = 1-eps/2');
%! assert(abs(elliptic3(-2.70143, 1 - eps/2, 0.9723) - (-1249.3300419938347)) < 1e-14 * 1249, 'elliptic3 reflection at m = 1-eps/2 (complete integral must be exact)');
%! assert(abs(carlsonRJ(0.1, 0.2, 1, 3.0) - 1.1311524759367163) < 5e-15 * 1.13, 'RJ series E3 coefficient');
%! assert(abs(carlsonRJ(0.292, 0.646, 1, 1.354) - 1.2806109121365949) < 5e-15 * 1.28, 'RJ series E3 coefficient');
%! [~, ~, J] = ellipticBDJ(1, 0.5, 1);
%! assert(abs(J - 0.64877476917835824) < 5e-15, 'J at n = 1 (was NaN)');
%! [~, ~, J] = ellipticBDJ(0.5, 0.5, 1.5);
%! assert(abs(J - 0.052791966372572887) < 5e-16, 'J at n = 1.5 below the pole');
%! err = '';
%! try, ellipticBDJ(1, 0.5, 1.5); catch e, err = e.message; end
%! assert(~isempty(strfind(err, 'principal')), 'J at n = 1.5 beyond the pole must error, not return complex');
%! assert(abs(elliptic3(1, 0.5, -0.5) - 0.9560406633267465) < 5e-16, 'elliptic3 with c = -0.5');
%! assert(abs(elliptic3(1, 0.5, -3.0) - 0.66684868942035313) < 5e-16, 'elliptic3 with c = -3');
%! assert(abs(elliptic3(1, 0.5, -100.0) - 0.1523863772236308) < 5e-16, 'elliptic3 with c = -100 (Carlson branch)');
%! assert(abs(elliptic3(4, 0.9, -100.0) - 0.4921742710224714) < 5e-15, 'elliptic3 with c = -100, reduced phase');

%% ---------------------------------------------------------------------
%% V. Bulirsch cel (round 6b).  The old route through m = 1 - kc^2 lost kc
%%    entirely below ~1e-8 (cel1(1e-9) was Inf here, 2e6 in Python; the value
%%    is ln(4/kc) = 22.1), rejected kc > 1 (m < 0) and returned Inf for p < 0.
%%    Bulirsch's own algorithm is kc-native; p < 0 is the Cauchy principal
%%    value, equal to Re Pi(1-p | 1-kc^2) (mpmath).
%% ---------------------------------------------------------------------
%!test
%! clear
%! assert(abs(cel1(1e-9) - (22.109560198066302)) < 1e-15 * 22, 'cel1(1e-9) = ln(4/kc)+...');
%! assert(abs(cel1(1e-300) - (692.1618222593336)) < 1e-15 * 692, 'cel1(1e-300)');
%! assert(abs(cel1(2) - (1.0782578237498216)) < 1e-15, 'cel1(2): kc > 1 is m = -3');
%! assert(abs(cel(0.5, -0.5, 1, 1) - (-1.0782578237498216)) < 1e-15 * 1.1, 'cel with p < 0 = principal value Re Pi(1.5|0.75)');
%! assert(abs(cel(0.7, -5, 1, 1) - (-0.092277884964284496)) < 1e-15, 'cel with p = -5');
%! assert(abs(cel(1e-9, 0.3, 1.5, -0.7) - (-46.045402351061091)) < 1e-15 * 47, 'general cel at kc = 1e-9');
%! assert(cel(0.3, 1, 1, 1) == cel1(0.3) && cel(-0.3, 1, 1, 1) == cel1(0.3), 'cel depends on kc^2 only');
%! assert(isinf(cel1(0)) && cel1(0) > 0, 'cel1(0) = K(1) = Inf');
%! kc = [1e-12 0.3 0.9 2.5];  p = [-0.4 0.7 -3 1e-3];
%! assert(all(abs(cel(kc, p, 1, 0) + p .* cel(kc, p, 0, 1) - cel1(kc)) < 1e-14 .* max(1, cel1(kc))), 'cel(1,0) + p cel(0,1) = K');

%% ---------------------------------------------------------------------
%% W. Round 6c (sweep of the remaining outputs).  mpmath at exact doubles.
%%    - ellipticBDJ formed Delta^2 = 1 - m sin^2, which cancels near phi = pi/2
%%      as m -> 1; jacobiEDJ also took am(u) at |u| ~ 1e3 before reducing, so
%%      D_u(1520|1-1e-8) was off by 3e-10.  Both are at the eps*|u| floor now.
%%    - arclength_ellipse branched with if(a<b) on arrays (all-elements
%%      semantics): mixed arrays fell through to the circle formula.
%% ---------------------------------------------------------------------
%!test
%! clear
%! [~, Dv] = ellipticBDJ(pi/2 - 2e-4, 1 - 1e-8);
%! assert(abs(Dv - (8.1529993379146678)) < 2e-12 * 9, 'D(pi/2-2e-4 | 1-1e-8): Delta^2 must not cancel');
%! [Eu, Du] = jacobiEDJ(1520.3427441800743, 0.999999990458008);
%! assert(abs(Eu - 143.00000694616405) < 2e-12 && abs(Du - 1377.3427503765037) < 2e-12, 'jacobiEDJ at u = 1520 near m = 1 (was 3e-10 off)');
%! v = arclength_ellipse([5 785.9 3], [10 495.8 3], [0 5.279 0], [1 -6.134 2]);
%! assert(all(abs(v - [8.8662512353670695 -7494.1448816975323 6]) < 1e-14 * [9 7495 6]), 'arclength_ellipse elementwise on mixed arrays');
%! assert(abs(arclength_ellipse(5, 10, [0 0.5], [1 1.5])(1) - 8.8662512353670695) < 1e-14, 'scalar axes with array angles');

%% ---------------------------------------------------------------------
%% X. inversenomeq / nome2m near q = 1.  Above q_max = 0.7789534 the true
%%    1-m = m(exp(pi^2/ln q)) is below eps/2, so the correctly rounded double
%%    is exactly 1; the unconverged 30-term series returned m > 1 (1.034 at
%%    q = 0.999, 1+1.6e-15 at q = 0.9).  nome2m captured its whole input array
%%    in the fzero objective and errored on any array.
%% ---------------------------------------------------------------------
%!test
%! clear
%! assert(isequal(inversenomeq([0.78 0.9 0.999]), [1 1 1]), 'inversenomeq rounds to exactly 1 above q_max');
%! assert(all(inversenomeq([0.5 0.7 0.7789]) <= 1), 'no overshoot below q_max');
%! assert(abs(inversenomeq(0.5) - 0.99998952213731039) < 4e-16, 'inversenomeq(0.5) (mpmath mfrom)');
%! q = [1e-12 0.05 0.3 0.5];          % q > 0.5 puts 1-m below 1e-5, where q(m) is ill-conditioned in double
%! m = nome2m(q);
%! assert(isequal(size(m), size(q)) && max(abs(nomeq(m) - q) ./ q) < 1e-11, 'nome2m on an array, round trip');
%! assert(nome2m(0.999) == 1, 'nome2m near q = 1');

%% ---------------------------------------------------------------------
%% Y. elliptic12i period term just below pi/2, and elliptic123 for m > 1.
%%    lambda = (-1)^k lambda + pi*ceil(phi/pi - 0.5 + eps) counted the period
%%    from a separately rounded quantity: for phi within a few ulps below
%%    pi/2 (asin(sqrt(3)) is 2 ulps below) Re F came out 3K instead of K,
%%    which elliptic123 inherited as K(3) = 3.003 (mpmath: 1.001).  The
%%    complete m > 1 case now uses the DLMF 19.7.3 closed forms instead of
%%    evaluating elliptic12i exactly on its branch point (sqrt(eps)-conditioned).
%% ---------------------------------------------------------------------
%!test
%! clear
%! F = elliptic12i(1.5707963267948961 + 0.5i, 1/3);
%! assert(abs(F - (1.7339168852579344 + 0.62666316872107993i)) < 2e-15 * 2, 'F two ulps below pi/2 (was 3K + ...)');
%! F = elliptic12i(asin(sqrt(3)), 1/3);
%! assert(abs(real(F) - 1.73391688525794) < 1e-7 && abs(imag(F) + 2.02895910274881) < 1e-7, 'F at the branch point (sqrt(eps)-conditioned there)');
%! [K, E] = elliptic123(3);
%! assert(abs(K - (1.0010773804561062 - 1.1714200841467699i)) < 1e-15 * 2 && abs(E - (0.47522393535101711 + 1.0130180585994313i)) < 1e-15 * 2, 'K(3), E(3) closed forms (were 3x off in the real part)');
%! [K, E] = elliptic123(pi/2, 3);
%! assert(abs(K - (1.0010773804561062 - 1.1714200841467699i)) < 1e-15 * 2 && abs(E - (0.47522393535101711 + 1.0130180585994313i)) < 1e-15 * 2, 'elliptic123(pi/2, 3) routes to the complete closed form');
%! [K, E] = elliptic123(5);
%! assert(abs(K - (0.74220623671119323 - 1.0094529099892116i)) < 3e-16 && abs(E - (0.36075866393790281 + 1.6257306716064185i)) < 3e-15, 'K(5), E(5)');
%! [F, E] = elliptic123(1.2, 3);
%! assert(abs(F - (1.0010773804561062 - 0.89956974520736591i)) < 1e-14 && abs(E - (0.47522393535101711 + 0.50673122232331459i)) < 1e-14, 'F(1.2|3), E(1.2|3) (m > 1, real part was 3x off)');
