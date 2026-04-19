%Test parallel vs serial alignment for all parallelized functions.
%Tests use a small chunk_size to force chunking even on small inputs.

%!test
%! % elliptic12: parallel output must match serial output
%! clear
%! [phi, alpha] = meshgrid(linspace(0, pi/2, 50), linspace(0, pi/2, 50));
%! u = phi(:).';
%! m = sin(alpha(:).').^2;
%! elliptic_config('parallel', false);
%! [F_s, E_s, Z_s] = elliptic12(u, m);
%! elliptic_config('parallel', true);
%! elliptic_config('chunk_size', 100);
%! [F_p, E_p, Z_p] = elliptic12(u, m);
%! elliptic_config('parallel', false);
%! assert(max(abs(F_p - F_s)) < 1e-14, 'elliptic12 F: parallel/serial mismatch');
%! assert(max(abs(E_p - E_s)) < 1e-14, 'elliptic12 E: parallel/serial mismatch');
%! assert(max(abs(Z_p - Z_s)) < 1e-14, 'elliptic12 Z: parallel/serial mismatch');

%!test
%! % elliptic3: parallel output must match serial output
%! clear
%! [phi, alpha, cv] = meshgrid(linspace(0, pi/2, 15), linspace(0, pi/2, 15), linspace(0, 0.9, 5));
%! u = phi(:).';
%! m = sin(alpha(:).').^2;
%! c = cv(:).';
%! elliptic_config('parallel', false);
%! Pi_s = elliptic3(u, m, c);
%! elliptic_config('parallel', true);
%! elliptic_config('chunk_size', 100);
%! Pi_p = elliptic3(u, m, c);
%! elliptic_config('parallel', false);
%! assert(max(abs(Pi_p - Pi_s)) < 1e-14, 'elliptic3: parallel/serial mismatch');

%!test
%! % ellipj: parallel output must match serial output
%! clear
%! [phi, alpha] = meshgrid(linspace(0, 10, 50), linspace(0, pi/2, 50));
%! u = phi(:).';
%! m = sin(alpha(:).').^2;
%! elliptic_config('parallel', false);
%! [Sn_s, Cn_s, Dn_s, Am_s] = ellipj(u, m);
%! elliptic_config('parallel', true);
%! elliptic_config('chunk_size', 100);
%! [Sn_p, Cn_p, Dn_p, Am_p] = ellipj(u, m);
%! elliptic_config('parallel', false);
%! assert(max(abs(Sn_p - Sn_s)) < 1e-14, 'ellipj Sn: parallel/serial mismatch');
%! assert(max(abs(Cn_p - Cn_s)) < 1e-14, 'ellipj Cn: parallel/serial mismatch');
%! assert(max(abs(Dn_p - Dn_s)) < 1e-14, 'ellipj Dn: parallel/serial mismatch');
%! assert(max(abs(Am_p - Am_s)) < 1e-14, 'ellipj Am: parallel/serial mismatch');

%!test
%! % jacobiThetaEta: parallel output must match serial output
%! clear
%! [phi, qv] = meshgrid(linspace(0, pi, 50), linspace(0.01, 0.9, 50));
%! u = phi(:).';
%! m = qv(:).'.^2;  % use m in (0,1) range
%! elliptic_config('parallel', false);
%! [Th_s, H_s] = jacobiThetaEta(u, m);
%! elliptic_config('parallel', true);
%! elliptic_config('chunk_size', 100);
%! [Th_p, H_p] = jacobiThetaEta(u, m);
%! elliptic_config('parallel', false);
%! assert(max(abs(Th_p - Th_s)) < 1e-14, 'jacobiThetaEta Th: parallel/serial mismatch');
%! assert(max(abs(H_p  - H_s )) < 1e-14, 'jacobiThetaEta H: parallel/serial mismatch');

%!test
%! % Verify serial fallback when parallel is disabled (no worker overhead)
%! clear
%! u = linspace(0, pi/2, 200);
%! m = linspace(0, 0.99, 200);
%! elliptic_config('parallel', false);
%! n = get_nworkers();
%! assert(n == 0, 'get_nworkers must return 0 when parallel is disabled');
%! [F1, E1] = elliptic12(u, m);
%! [F2, E2] = elliptic12(u, m);
%! assert(isequal(F1, F2), 'Serial results must be deterministic');
%! assert(isequal(E1, E2), 'Serial results must be deterministic');
