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

% ---------------------------------------------------------------------
% Chunking path exercised WITHOUT the parallel package: a temporary dir
% shadows get_nworkers (3 workers) and provides a serial parcellfun with
% the same calling convention.  Covers N below, at and above exact
% multiples of chunk_size -- the exact-multiple case re-entered the
% dispatch from inside par_worker without bound and crashed Octave.
% ---------------------------------------------------------------------
%!test
%! d = tempname(); mkdir(d);
%! fid = fopen(fullfile(d, 'get_nworkers.m'), 'w');
%! fprintf(fid, 'function n = get_nworkers()\nif ~elliptic_config(''parallel''), n = 0; else n = 3; end\n'); fclose(fid);
%! fid = fopen(fullfile(d, 'parcellfun.m'), 'w');
%! fprintf(fid, ['function varargout = parcellfun(nproc, fn, varargin)\n' ...
%!     'uo = true; args = varargin; k = find(strcmp(args, ''UniformOutput''), 1);\n' ...
%!     'if ~isempty(k), uo = args{k+1}; args(k:k+1) = []; end\n' ...
%!     'n = numel(args{1}); out = cell(1, n);\n' ...
%!     'for i = 1:n, a = cell(1, numel(args)); for j = 1:numel(args), a{j} = args{j}{i}; end; out{i} = fn(a{:}); end\n' ...
%!     'if uo, varargout{1} = [out{:}]; else varargout{1} = out; end\n']); fclose(fid);
%! addpath(d);
%! old_par = elliptic_config('parallel'); old_cs = elliptic_config('chunk_size');
%! unwind_protect
%!     cs = 500;
%!     for N = [cs-1, cs, 2*cs, 3*cs, 3*cs+7]
%!         rand('seed', N); u = rand(1,N)*40-20; m = rand(1,N)*0.98+0.01; c = rand(1,N)*0.9;
%!         elliptic_config('parallel', false);
%!         [F1,E1] = elliptic12(u,m); s1 = ellipj(u,m); P1 = elliptic3(u,m,c); [B1,D1,J1] = ellipticBDJ(u,m,c);
%!         elliptic_config('parallel', true); elliptic_config('chunk_size', cs);
%!         [F2,E2] = elliptic12(u,m); s2 = ellipj(u,m); P2 = elliptic3(u,m,c); [B2,D2,J2] = ellipticBDJ(u,m,c);
%!         assert(isequal(size(F2), size(F1)) && max(abs(F2-F1)) == 0 && max(abs(E2-E1)) == 0, 'elliptic12 chunked != serial at N=%d', N);
%!         assert(max(abs(s2-s1)) == 0, 'ellipj chunked != serial at N=%d', N);
%!         assert(max(abs(P2-P1)) == 0, 'elliptic3 chunked != serial at N=%d', N);
%!         assert(max(abs(J2-J1)) == 0 && max(abs(B2-B1)) == 0, 'ellipticBDJ chunked != serial at N=%d', N);
%!     end
%!     [ph, al] = meshgrid(linspace(0.1, 3, 40), linspace(0.05, 0.9, 30));   % matrix through the chunked path
%!     Fm = elliptic12(ph, al); elliptic_config('parallel', false); Fs = elliptic12(ph, al);
%!     assert(isequal(size(Fm), [30 40]) && max(abs(Fm(:)-Fs(:))) == 0, 'matrix shape/values through chunking');
%!     assert(elliptic_config('parallel') == false, 'par_worker must restore the parallel flag');
%! unwind_protect_cleanup
%!     elliptic_config('parallel', old_par); elliptic_config('chunk_size', old_cs);
%!     rmpath(d); confirm_recursive_rmdir(false, 'local'); rmdir(d, 's');
%! end_unwind_protect
