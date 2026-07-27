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
