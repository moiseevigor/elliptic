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

% ------------------------------------------------------------------
%  Group 4: Hardcoded reference tables
%
%  ellipticBDJ: B(phi|m), D(phi|m), J(phi,n=0.3|m) on a 3-row x 4-col
%  meshgrid: m in {0.2, 0.5, 0.8}, phi in {pi/6, pi/4, pi/3, pi/2}.
%  Cross-checked: max|B+D - F| < 1e-15, max|B+(1-m)*D - E| < 2e-14.
%
%  jacobiEDJ: E_u, D_u, J_u(n=0.3) at u = {K/4, K/2, 3K/4, K} for
%  m in {0.2, 0.5, 0.8}.  Cross-checked: max|E_u-(u-m*D_u)| < 1e-15.
% ------------------------------------------------------------------

%!test
%! % ellipticBDJ B reference table: rows = m, cols = phi
%! phi_v = [pi/6, pi/4, pi/3, pi/2];
%! m_v   = [0.2;  0.5;  0.8 ];
%! [PHI, M] = meshgrid(phi_v, m_v);
%! expectedB = [ ...
%!   0.4822319097177440, 0.6529667634636220, 0.7569756041614870, 0.8066808960371530; ...
%!   0.4884759118954500, 0.6703551321060780, 0.7874738572402080, 0.8472130847939790; ...
%!   0.4952048001091030, 0.6909260550736680, 0.8281306941693560, 0.9088110737045840];
%! [B, ~] = ellipticBDJ(PHI, M);
%! assert(norm(B - expectedB) < 1e-12, 'ellipticBDJ B reference table failed');

%!test
%! % ellipticBDJ D reference table: rows = m, cols = phi
%! phi_v = [pi/6, pi/4, pi/3, pi/2];
%! m_v   = [0.2;  0.5;  0.8 ];
%! [PHI, M] = meshgrid(phi_v, m_v);
%! expectedD = [ ...
%!   4.600298959653437e-02, 1.474128859214848e-01, 3.234022020908615e-01, 8.529427025733753e-01; ...
%!   4.714682090995279e-02, 1.556627441431676e-01, 3.549552008055692e-01, 1.006861592507392e+00; ...
%!   4.839998016351204e-02, 1.658742182004513e-01, 4.016987480555824e-01, 1.348394253116269e+00];
%! [~, D] = ellipticBDJ(PHI, M);
%! assert(norm(D - expectedD) < 1e-12, 'ellipticBDJ D reference table failed');

%!test
%! % ellipticBDJ J reference table (n=0.3): rows = m, cols = phi
%! phi_v = [pi/6, pi/4, pi/3, pi/2];
%! m_v   = [0.2;  0.5;  0.8 ];
%! n     = 0.3;
%! [PHI, M] = meshgrid(phi_v, m_v);
%! expectedJ = [ ...
%!   4.823487054607033e-02, 1.630222842878820e-01, 3.807339316051869e-01, 1.112925198627194; ...
%!   4.944500507300179e-02, 1.723183252575266e-01, 4.189753046197145e-01, 1.321007148808583; ...
%!   5.077140114947837e-02, 1.838496301039670e-01, 4.759688499155832e-01, 1.788630886401574];
%! [~, ~, J] = ellipticBDJ(PHI, M, n .* ones(size(PHI)));
%! assert(norm(J - expectedJ) < 1e-12, 'ellipticBDJ J reference table failed');

%!test
%! % jacobiEDJ E_u reference table: rows = m, cols = u/K fraction
%! % u_frac in {1/4, 1/2, 3/4, 1} of K(m)
%! m_v = [0.2, 0.5, 0.8];
%! expectedEu = [ ...
%!   0.4103348919732966, 0.7973039335479685, 1.153380434524441, 1.489035058095853; ...
%!   0.4479288836359817, 0.8217685499305638, 1.110593848820260, 1.350643881047676; ...
%!   0.5213087503809286, 0.8656381644139404, 1.055687219771291, 1.178489924327838];
%! for mi = 1:3
%!   [K, ~] = ellipke(m_v(mi));
%!   u_v = [0.25*K, 0.5*K, 0.75*K, K];
%!   [Eu, ~] = jacobiEDJ(u_v, m_v(mi));
%!   assert(max(abs(Eu - expectedEu(mi,:))) < 1e-13, ...
%!       sprintf('jacobiEDJ Eu table failed at m=%.1f', m_v(mi)));
%! end

%!test
%! % jacobiEDJ D_u reference table: rows = m, cols = u/K fraction
%! m_v = [0.2, 0.5, 0.8];
%! expectedDu = [ ...
%!   0.02285503839667675, 0.1625393287864774, 0.4566863221672750, 0.8529427025733753; ...
%!   0.03117957137872241, 0.2105375774402440, 0.5599243183115373, 1.006861592507392;  ...
%!   0.05374072665535606, 0.3287056237456081, 0.7965209691804365, 1.348394253116269];
%! for mi = 1:3
%!   [K, ~] = ellipke(m_v(mi));
%!   u_v = [0.25*K, 0.5*K, 0.75*K, K];
%!   [~, Du] = jacobiEDJ(u_v, m_v(mi));
%!   assert(max(abs(Du - expectedDu(mi,:))) < 1e-13, ...
%!       sprintf('jacobiEDJ Du table failed at m=%.1f', m_v(mi)));
%! end

%!test
%! % jacobiEDJ J_u reference table (n=0.3): rows = m, cols = u/K fraction
%! m_v = [0.2, 0.5, 0.8];
%! expectedJu = [ ...
%!   0.02354858117828542, 0.1809070499958038, 0.5570048420102275, 1.112925198627194; ...
%!   0.03233444406148652, 0.2379084919244906, 0.6917734853508237, 1.321007148808583; ...
%!   0.05656046648375568, 0.3828083834319858, 1.007421432719784,  1.788630886401574];
%! for mi = 1:3
%!   [K, ~] = ellipke(m_v(mi));
%!   u_v = [0.25*K, 0.5*K, 0.75*K, K];
%!   [~, ~, Ju] = jacobiEDJ(u_v, m_v(mi), 0.3 .* ones(size(u_v)));
%!   assert(max(abs(Ju - expectedJu(mi,:))) < 1e-13, ...
%!       sprintf('jacobiEDJ Ju table failed at m=%.1f', m_v(mi)));
%! end
