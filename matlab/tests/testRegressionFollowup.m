% Regression coverage from the post-0d09740 deep audit.

%!test
%! % m=1: period reduction must not hide the first F pole or under-count E.
%! phi = [2, pi, 4, 10];
%! [F, E, Z] = elliptic12(phi, ones(size(phi)));
%! turns = floor((abs(phi) + pi/2) ./ pi);
%! expectedE = (-1).^turns .* sin(abs(phi)) + 2.*turns;
%! assert(all(isinf(F) & F > 0), 'F(phi|1) must diverge after crossing pi/2');
%! assert(max(abs(E-expectedE)) < 1e-14, 'E(phi|1) period accounting failed');
%! [Fn, En, Zn] = elliptic12(-phi, ones(size(phi)));
%! assert(all(isinf(Fn) & Fn < 0), 'negative F(phi|1) pole sign failed');
%! assert(max(abs(En+E)) < 1e-14 && max(abs(Zn+Z)) < 1e-14, 'm=1 parity failed');

%!test
%! % B, D and especially S retain their analytic limits for tiny m.
%! m = [0, 1e-20, 1e-16, 1e-12, 1e-8];
%! [B, D, S] = ellipticBD(m);
%! assert(max(abs(B-pi/4)) < 2e-8, 'B lost its m->0 limit');
%! assert(max(abs(D-pi/4)) < 2e-8, 'D lost its m->0 limit');
%! assert(max(abs(S-pi/16)) < 2e-8, 'S suffered small-m cancellation');

%!test
%! % Near-pole third-kind value anchored to scipy Carlson RF/RJ.
%! got = elliptic3(pi/2, 0.9, 0.999);
%! expected = 149.26048203240563;
%! assert(abs(got-expected) < 2e-12, ...
%!     'elliptic3 near-pole error: got %.17g expected %.17g', got, expected);

%!test
%! % Large u and the closest representable m<1: mpmath 50-digit anchor.
%! u = 1000000.123;
%! m = 1 - eps/2;
%! [sn, cn, dn] = ellipj(u, m);
%! assert(abs(sn-0.9999999999999987) < 2e-14, 'large-u sn lost phase');
%! assert(abs(cn-5.0691640447381745e-08) < 2e-14, 'large-u cn lost phase');
%! assert(abs(dn-5.177513605688685e-08) < 2e-14, 'large-u dn cancellation');

%!test
%! % Distinct m values must never be tolerance-grouped into one data point.
%! u = 1.56;
%! m = [0.999, 0.999+5e-12, 0.5, 0.5+5e-12];
%! [Fv, Ev] = elliptic12(u*ones(size(m)), m);
%! for k = 1:numel(m)
%!     [Fs, Es] = elliptic12(u, m(k));
%!     assert(abs(Fv(k)-Fs) < 2e-13 && abs(Ev(k)-Es) < 2e-13, ...
%!         'array evaluation substituted m at index %d', k);
%! end

%!test
%! % elliptic123 must route through the repaired public implementation.
%! phi = [0.4, 1.2, 2.0, 4.0];
%! m = 0.4*ones(size(phi));
%! [F123, E123] = elliptic123(phi, m);
%! [F, E] = elliptic12(phi, m);
%! assert(max(abs(F123-F)) < 2e-12 && max(abs(E123-E)) < 2e-12, ...
%!     'elliptic123 retained a stale private elliptic12 implementation');
