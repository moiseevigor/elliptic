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
