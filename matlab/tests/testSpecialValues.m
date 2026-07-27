% Tests for special / boundary input values not covered by the existing suite.
%
% New coverage:
%   elliptic12  – m=1 closed forms: F=atanh(sin φ), E=sin φ, Z=sin φ
%                 φ=π periodicity: F(π|m)=2K(m), E(π|m)=2E(m)
%                 φ<0 odd symmetry at boundary values m=0 and m=1
%   ellipj      – quarter-period K: sn=1, cn=0, dn=√(1−m), am=π/2
%                 half-period 2K: sn=0, cn=−1, dn=1, am=π
%                 full-period 4K: sn=0, cn=1, dn=1
%                 parity: sn odd, cn and dn even
%   elliptic3   – n=0 degenerate: Π(φ|m,0)=F(φ|m)
%                 φ↦−φ antisymmetry
%                 φ=π periodicity: Π(π|m,n)=2Π(π/2|m,n)
%   ellipticBD  – m→1 limit: B(m)→1, D(m)→∞
%   ellipticBDJ – φ=π periodicity: B(π|m)=2B(m), D(π|m)=2D(m)
%   carlsonRC   – RC(0,y)=π/(2√y) for untested y∈{9, 1/9}
%                 RC(1,1)=1 as isolated scalar check
%   carlsonRF   – RF(0,1,1)=π/2 (=K(0)); RF(1,1,1)=1
%   carlsonRD   – RD(1,1,1)=1
%   carlsonRJ   – RJ(1,1,1,1)=1
%   agm         – agm(a,0,a)=0 for a>0; agm(1,1,0)→1; agm(0,0,0)→0
%   nomeq       – nomeq(1)→1 boundary
%   cel1        – cel1(1)=π/2 (=K(0))
%   theta       – θ₂(π/2,m)=0 (all cos terms vanish)
%                 θ₁(π,m)=0  (pseudo-period: θ₁(z+π)=−θ₁(z))
%                 Jacobi derivative product: θ₁′(0)=θ₂(0)·θ₃(0)·θ₄(0)
%
% Run with:  test testSpecialValues

% =========================================================================
%  elliptic12 — m=1 closed forms
% =========================================================================

%!test
%! % F(φ|1) = atanh(sin φ)  for |φ| < π/2
%! phis = [pi/6, pi/4, pi/3];
%! for k = 1:3
%!     [F, ~, ~] = elliptic12(phis(k), 1);
%!     assert(abs(F - atanh(sin(phis(k)))) < 1e-12, ...
%!         sprintf('F(%.4f|m=1): got %g, expected %g', phis(k), F, atanh(sin(phis(k)))));
%! end

%!test
%! % E(φ|1) = sin φ  for φ ∈ [0, π/2]
%! phis = [0, pi/6, pi/4, pi/3, pi/2];
%! for k = 1:5
%!     [~, E, ~] = elliptic12(phis(k), 1);
%!     assert(abs(E - sin(phis(k))) < 1e-12, ...
%!         sprintf('E(%.4f|m=1): got %g, expected %g', phis(k), E, sin(phis(k))));
%! end

%!test
%! % Z(φ|1) = sin φ  [Z = E − (E(m)/K(m))·F; E(1)/K(1)→0 so Z→E=sin φ]
%! phis = [pi/6, pi/4, pi/3];
%! for k = 1:3
%!     [~, ~, Z] = elliptic12(phis(k), 1);
%!     assert(abs(Z - sin(phis(k))) < 1e-12, ...
%!         sprintf('Z(%.4f|m=1): got %g, expected %g', phis(k), Z, sin(phis(k))));
%! end

% =========================================================================
%  elliptic12 — φ=π periodicity
% =========================================================================

%!test
%! % F(π|m) = 2 K(m)
%! m_vals = [0.2, 0.5, 0.7, 0.9];
%! for k = 1:4
%!     m = m_vals(k);
%!     [K, ~]    = ellipke(m);
%!     [F_pi, ~] = elliptic12(pi, m);
%!     assert(abs(F_pi - 2*K) < 1e-12, ...
%!         sprintf('F(π|m=%.1f): got %g, expected %g', m, F_pi, 2*K));
%! end

%!test
%! % E(π|m) = 2 E(m)
%! m_vals = [0.2, 0.5, 0.7, 0.9];
%! for k = 1:4
%!     m = m_vals(k);
%!     [~, E_c]  = ellipke(m);
%!     [~, E_pi] = elliptic12(pi, m);
%!     assert(abs(E_pi - 2*E_c) < 1e-12, ...
%!         sprintf('E(π|m=%.1f): got %g, expected %g', m, E_pi, 2*E_c));
%! end

% =========================================================================
%  elliptic12 — negative φ odd symmetry at boundary m values
% =========================================================================

%!test
%! % F(−φ|m)=−F(φ|m) and E(−φ|m)=−E(φ|m) at m=0 and m=1
%! phi = pi/5;
%! for m = [0, 1]
%!     [F_p, E_p] = elliptic12( phi, m);
%!     [F_n, E_n] = elliptic12(-phi, m);
%!     assert(abs(F_n + F_p) < 1e-12, sprintf('F odd symmetry failed at m=%g', m));
%!     assert(abs(E_n + E_p) < 1e-12, sprintf('E odd symmetry failed at m=%g', m));
%! end

% =========================================================================
%  ellipj — quarter-period values at u = K(m)
% =========================================================================

%!test
%! % sn(K|m)=1, cn(K|m)=0, dn(K|m)=√(1−m), am(K|m)=π/2
%! m_vals = [0.1, 0.3, 0.5, 0.7, 0.9];
%! for k = 1:5
%!     m = m_vals(k);
%!     [K, ~] = ellipke(m);
%!     [Sn, Cn, Dn, Am] = ellipj(K, m);
%!     assert(abs(Sn - 1.0)         < 1e-12, sprintf('sn(K|m=%.1f) != 1',       m));
%!     assert(abs(Cn - 0.0)         < 1e-12, sprintf('cn(K|m=%.1f) != 0',       m));
%!     assert(abs(Dn - sqrt(1-m))   < 1e-12, sprintf('dn(K|m=%.1f) != sqrt(mc)',m));
%!     assert(abs(Am - pi/2)        < 1e-12, sprintf('am(K|m=%.1f) != pi/2',    m));
%! end

% =========================================================================
%  ellipj — half-period values at u = 2K(m)
% =========================================================================

%!test
%! % sn(2K|m)=0, cn(2K|m)=−1, dn(2K|m)=1, am(2K|m)=π
%! m_vals = [0.1, 0.3, 0.5, 0.7, 0.9];
%! for k = 1:5
%!     m = m_vals(k);
%!     [K, ~] = ellipke(m);
%!     [Sn, Cn, Dn, Am] = ellipj(2*K, m);
%!     assert(abs(Sn - 0.0)  < 1e-10, sprintf('sn(2K|m=%.1f) != 0',  m));
%!     assert(abs(Cn + 1.0)  < 1e-10, sprintf('cn(2K|m=%.1f) != -1', m));
%!     assert(abs(Dn - 1.0)  < 1e-10, sprintf('dn(2K|m=%.1f) != 1',  m));
%!     assert(abs(Am - pi)   < 1e-10, sprintf('am(2K|m=%.1f) != pi', m));
%! end

% =========================================================================
%  ellipj — full-period restoration at u = 4K(m)
% =========================================================================

%!test
%! % sn(4K|m)=0, cn(4K|m)=1, dn(4K|m)=1  (tolerance relaxed for 4K accumulation)
%! m_vals = [0.3, 0.5, 0.7];
%! for k = 1:3
%!     m = m_vals(k);
%!     [K, ~] = ellipke(m);
%!     [Sn, Cn, Dn, ~] = ellipj(4*K, m);
%!     assert(abs(Sn - 0.0) < 1e-9, sprintf('sn(4K|m=%.1f) != 0', m));
%!     assert(abs(Cn - 1.0) < 1e-9, sprintf('cn(4K|m=%.1f) != 1', m));
%!     assert(abs(Dn - 1.0) < 1e-9, sprintf('dn(4K|m=%.1f) != 1', m));
%! end

% =========================================================================
%  ellipj — parity: sn odd, cn and dn even
% =========================================================================

%!test
%! % sn(−u|m)=−sn(u|m),  cn(−u|m)=cn(u|m),  dn(−u|m)=dn(u|m)
%! u_vals = [0.5, 1.2, 2.0];
%! m_vals = [0.0, 0.3, 0.7, 1.0];
%! for mi = 1:4
%!     m = m_vals(mi);
%!     for ui = 1:3
%!         u = u_vals(ui);
%!         [Sn_p, Cn_p, Dn_p] = ellipj( u, m);
%!         [Sn_n, Cn_n, Dn_n] = ellipj(-u, m);
%!         assert(abs(Sn_n + Sn_p) < 1e-14, sprintf('sn not odd  u=%.1f m=%.1f', u, m));
%!         assert(abs(Cn_n - Cn_p) < 1e-14, sprintf('cn not even u=%.1f m=%.1f', u, m));
%!         assert(abs(Dn_n - Dn_p) < 1e-14, sprintf('dn not even u=%.1f m=%.1f', u, m));
%!     end
%! end

% =========================================================================
%  elliptic3 — n=0 degenerate case: Π(φ|m,0) = F(φ|m)
% =========================================================================

%!test
%! phi_vals = [pi/6, pi/4, pi/3, pi/2];
%! m_vals   = [0.2, 0.5, 0.8];
%! for mi = 1:3
%!     for pi_i = 1:4
%!         m   = m_vals(mi);
%!         phi = phi_vals(pi_i);
%!         [F, ~] = elliptic12(phi, m);
%!         Pi0    = elliptic3(phi, m, 0);
%!         assert(abs(Pi0 - F) < 1e-11, ...
%!             sprintf('Pi(%.3f|m=%.1f,n=0)=%g != F=%g', phi, m, Pi0, F));
%!     end
%! end

% =========================================================================
%  elliptic3 — φ↦−φ antisymmetry
% =========================================================================

%!test
%! % Π(−φ|m,n) = −Π(φ|m,n)
%! phi   = 0.7;
%! n_vals = [0.2, 0.5, 0.8];
%! m_vals = [0.3, 0.6];
%! for mi = 1:2
%!     for ni = 1:3
%!         m = m_vals(mi);
%!         n = n_vals(ni);
%!         Pi_p = elliptic3( phi, m, n);
%!         Pi_n = elliptic3(-phi, m, n);
%!         assert(abs(Pi_n + Pi_p) < 1e-12, ...
%!             sprintf('Pi not odd in phi  m=%.1f n=%.1f', m, n));
%!     end
%! end

% =========================================================================
%  elliptic3 — φ=π periodicity: Π(π|m,n) = 2 Π(π/2|m,n)
% =========================================================================

%!test
%! m_vals = [0.2, 0.5, 0.7];
%! n_vals = [0.2, 0.5];
%! for mi = 1:3
%!     for ni = 1:2
%!         m = m_vals(mi);
%!         n = n_vals(ni);
%!         Pi_half = elliptic3(pi/2, m, n);
%!         Pi_full = elliptic3(pi,   m, n);
%!         assert(abs(Pi_full - 2*Pi_half) < 1e-11, ...
%!             sprintf('Pi(pi)!=2Pi(pi/2)  m=%.1f n=%.1f', m, n));
%!     end
%! end

% =========================================================================
%  ellipticBD — m→1 limit: B(m)→1, D(m)→∞
% =========================================================================

%!test
%! % B(m) → 1 as m → 1  [B=(E-(1-m)K)/m, (1-m)K→0, E→1]
%! m_vals = [0.9, 0.99, 0.999, 0.9999];
%! [B, ~, ~] = ellipticBD(m_vals);
%! assert(abs(B(end) - 1.0) < 5e-4, ...
%!     sprintf('B(0.9999) = %g, expected ≈1', B(end)));
%! assert(all(diff(B) > 0), 'B(m) should be increasing toward m=1');

%!test
%! % D(m) = (K-E)/m grows monotonically as m → 1
%! m_vals = [0.9, 0.99, 0.999, 0.9999];
%! D_vals = zeros(1, 4);
%! for k = 1:4
%!     [~, D_vals(k), ~] = ellipticBD(m_vals(k));
%! end
%! assert(all(diff(D_vals) > 0), 'D(m) should be monotone increasing toward m=1');
%! assert(D_vals(end) > D_vals(1), ...
%!     sprintf('D(0.9999)=%g should exceed D(0.9)=%g', D_vals(end), D_vals(1)));

% =========================================================================
%  ellipticBDJ — φ=π periodicity: B(π|m)=2B(m), D(π|m)=2D(m)
% =========================================================================

%!test
%! % B(π|m) = 2 B_complete(m),  D(π|m) = 2 D_complete(m)
%! m_vals = [0.2, 0.5, 0.7, 0.9];
%! for k = 1:4
%!     m = m_vals(k);
%!     [B_c, D_c, ~] = ellipticBD(m);
%!     [B_pi, D_pi]  = ellipticBDJ(pi, m);
%!     assert(abs(B_pi - 2*B_c) < 1e-11, ...
%!         sprintf('B(pi|m=%.1f): got %g, expected %g', m, B_pi, 2*B_c));
%!     assert(abs(D_pi - 2*D_c) < 1e-11, ...
%!         sprintf('D(pi|m=%.1f): got %g, expected %g', m, D_pi, 2*D_c));
%! end

% =========================================================================
%  carlsonRC — additional RC(0,y)=π/(2√y) values
% =========================================================================

%!test
%! % RC(0, 9) = π/6
%! RC = carlsonRC(0, 9);
%! assert(abs(RC - pi/6) < 1e-14, sprintf('RC(0,9) = %g, expected pi/6 = %g', RC, pi/6));

%!test
%! % RC(0, 1/9) = 3π/2
%! RC = carlsonRC(0, 1/9);
%! assert(abs(RC - 3*pi/2) < 1e-14, sprintf('RC(0,1/9) = %g, expected 3pi/2 = %g', RC, 3*pi/2));

%!test
%! % RC(1, 1) = 1  (isolated scalar; the x=1 case from the RC(x,x)=1/√x array test)
%! RC = carlsonRC(1, 1);
%! assert(abs(RC - 1.0) < 1e-14, sprintf('RC(1,1) = %g, expected 1', RC));

% =========================================================================
%  carlsonRF — RF(0,1,1)=π/2 and RF(1,1,1)=1
% =========================================================================

%!test
%! % RF(0, 1, 1) = π/2  [RF(0,y,y)=RC(0,y)=π/(2√y); at y=1 → π/2 = K(0)]
%! RF = carlsonRF(0, 1, 1);
%! assert(abs(RF - pi/2) < 1e-13, sprintf('RF(0,1,1) = %g, expected pi/2', RF));

%!test
%! % RF(1, 1, 1) = 1  [RF(x,x,x) = x^(-1/2); isolated scalar at x=1]
%! RF = carlsonRF(1, 1, 1);
%! assert(abs(RF - 1.0) < 1e-13, sprintf('RF(1,1,1) = %g, expected 1', RF));

%!test
%! % RF(0, 4, 4) = π/4  [RC(0,4)=π/(2·2)=π/4]
%! RF = carlsonRF(0, 4, 4);
%! assert(abs(RF - pi/4) < 1e-13, sprintf('RF(0,4,4) = %g, expected pi/4', RF));

% =========================================================================
%  carlsonRD — RD(1,1,1)=1
% =========================================================================

%!test
%! % RD(1, 1, 1) = 1  [RD(x,x,x)=x^(-3/2); isolated scalar at x=1]
%! RD = carlsonRD(1, 1, 1);
%! assert(abs(RD - 1.0) < 1e-13, sprintf('RD(1,1,1) = %g, expected 1', RD));

% =========================================================================
%  carlsonRJ — RJ(1,1,1,1)=1
% =========================================================================

%!test
%! % RJ(1,1,1,1) = 1  [RJ(x,y,z,z)=RD(x,y,z); at x=y=z=p=1: RD(1,1,1)=1]
%! RJ = carlsonRJ(1, 1, 1, 1);
%! assert(abs(RJ - 1.0) < 1e-13, sprintf('RJ(1,1,1,1) = %g, expected 1', RJ));

% =========================================================================
%  AGM — agm(a,0,a)→0, agm(1,1,0)→1, agm(0,0,0)→0
% =========================================================================

%!test
%! % agm(a, 0, a) converges to 0  [AGM of a and 0 is 0 for any a>0]
%! a_vals = [0.5, 1, 2, pi];
%! for k = 1:4
%!     a = a_vals(k);
%!     [av, bv, ~, nv] = agm(a, 0, a);
%!     mn = double(max(nv));
%!     assert(abs(av(mn+1)) < 1e-14, ...
%!         sprintf('agm(%g,0) should converge to 0, got %g', a, av(mn+1)));
%! end

%!test
%! % agm(1, 1, 0) converges to 1  [AGM of equal args is that arg]
%! [av, bv, ~, nv] = agm(1, 1, 0);
%! mn = double(max(nv));
%! assert(abs(av(mn+1) - 1.0) < 1e-15, ...
%!     sprintf('agm(1,1) should be 1, got %g', av(mn+1)));

%!test
%! % agm(0, 0, 0) converges to 0
%! [av, bv, ~, nv] = agm(0, 0, 0);
%! mn = double(max(nv));
%! assert(abs(av(mn+1)) < 1e-15, 'agm(0,0) should be 0');

%!test
%! % AGM is symmetric: agm(a,b) = agm(b,a)
%! pairs = [2.0, 1.0; 0.3, 0.7; pi, 1.0];
%! for k = 1:3
%!     a = pairs(k,1); b = pairs(k,2);
%!     [av1, ~, ~, nv1] = agm(a, b, a-b);
%!     [av2, ~, ~, nv2] = agm(b, a, a-b);
%!     mn1 = double(max(nv1)); mn2 = double(max(nv2));
%!     assert(abs(av1(mn1+1) - av2(mn2+1)) < 1e-14, ...
%!         sprintf('AGM not symmetric: agm(%g,%g)=%g vs agm(%g,%g)=%g', ...
%!             a,b,av1(mn1+1), b,a,av2(mn2+1)));
%! end

% =========================================================================
%  nomeq — boundary values
% =========================================================================

%!test
%! % nomeq(1) should approach 1  [q=exp(-pi*K'/K) → 1 as K'/K → 0]
%! q = nomeq(1);
%! assert(abs(q - 1.0) < 1e-4, sprintf('nomeq(1) = %g, expected ≈1', q));

%!test
%! % nomeq is monotone: nomeq(0.9999) > nomeq(0.999) > ... > nomeq(0)
%! m_vals = [0, 0.1, 0.5, 0.9, 0.999, 0.9999];
%! q_vals = nomeq(m_vals);
%! assert(all(diff(q_vals) > 0), 'nomeq should be strictly increasing');

% =========================================================================
%  Bulirsch cel1 — kc=1 gives K(0)=π/2
% =========================================================================

%!test
%! % cel1(1) = K(1 - 1^2) = K(0) = π/2
%! val = cel1(1);
%! assert(abs(val - pi/2) < 1e-13, sprintf('cel1(1) = %g, expected pi/2', val));

% =========================================================================
%  theta — θ₂(π/2, m) = 0  [cos((2n+1)π/2)=0 for all n]
% =========================================================================

%!test
%! m_vals = [0.1, 0.3, 0.5, 0.7, 0.9];
%! for k = 1:5
%!     m = m_vals(k);
%!     [th2, ~] = theta_prime(2, pi/2, m);
%!     assert(abs(th2) < 1e-12, ...
%!         sprintf('theta_2(pi/2 | m=%.1f) = %g, expected 0', m, th2));
%! end

% =========================================================================
%  theta — θ₁(π, m) = 0  [θ₁(z+π)=−θ₁(z), so θ₁(π)=−θ₁(0)=0]
% =========================================================================

%!test
%! m_vals = [0.1, 0.3, 0.5, 0.7, 0.9];
%! for k = 1:5
%!     m = m_vals(k);
%!     [th1, ~] = theta_prime(1, pi, m);
%!     assert(abs(th1) < 1e-12, ...
%!         sprintf('theta_1(pi | m=%.1f) = %g, expected 0', m, th1));
%! end

% =========================================================================
%  theta — Jacobi derivative product: θ₁′(0) = θ₂(0)·θ₃(0)·θ₄(0)
% =========================================================================

%!test
%! m_vals = [0.1, 0.3, 0.5, 0.7, 0.9];
%! for k = 1:5
%!     m = m_vals(k);
%!     [~,  dth1] = theta_prime(1, 0, m);
%!     [th2, ~  ] = theta_prime(2, 0, m);
%!     [th3, ~  ] = theta_prime(3, 0, m);
%!     [th4, ~  ] = theta_prime(4, 0, m);
%!     product = th2 * th3 * th4;
%!     tol = 1e-11 * max(abs(product), 1);
%!     assert(abs(dth1 - product) < tol, ...
%!         sprintf('theta_1''(0) = %g, product = %g  (m=%.1f)', dth1, product, m));
%! end
