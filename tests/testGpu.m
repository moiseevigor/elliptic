%Test GPU acceleration path against serial CPU results.
%Each test skips gracefully when no GPU hardware is present.
%Enable with: elliptic_config('gpu', true)

%!function tf = gpu_available_for_test()
%! % Returns true if a GPU is physically present (ignores elliptic_config).
%! tf = false;
%! if exist('OCTAVE_VERSION', 'builtin')
%!     % Octave: try the 'ocl' Forge package (OpenCL)
%!     try
%!         pkg('load', 'ocl');
%!         tmp = gpuArray(1.0);
%!         gather(tmp);
%!         tf = true;
%!     catch
%!         tf = false;
%!     end
%!     return;
%! end
%! if exist('canUseGPU', 'builtin') || exist('canUseGPU', 'file')
%!     tf = canUseGPU();
%! elseif exist('gpuDeviceCount', 'file')
%!     tf = gpuDeviceCount() > 0;
%! end

%!test
%! % elliptic12: GPU output must match serial CPU output
%! clear
%! elliptic_config('gpu', false);
%! if ~gpu_available_for_test(), disp('SKIP: no GPU'); return; end
%! [phi, alpha] = meshgrid(linspace(0.01, pi/2, 40), linspace(0.01, pi/2, 40));
%! u = phi(:).'; m = sin(alpha(:).').^2;
%! [F_s, E_s, Z_s] = elliptic12(u, m);
%! elliptic_config('gpu', true);
%! [F_g, E_g, Z_g] = elliptic12(u, m);
%! elliptic_config('gpu', false);
%! assert(max(abs(F_g - F_s)) < 1e-12, 'elliptic12 F: GPU/serial mismatch');
%! assert(max(abs(E_g - E_s)) < 1e-12, 'elliptic12 E: GPU/serial mismatch');
%! assert(max(abs(Z_g - Z_s)) < 1e-12, 'elliptic12 Z: GPU/serial mismatch');

%!test
%! % elliptic3: GPU output must match serial CPU output
%! clear
%! elliptic_config('gpu', false);
%! if ~gpu_available_for_test(), disp('SKIP: no GPU'); return; end
%! [phi, alpha, cv] = meshgrid(linspace(0, pi/2, 20), linspace(0, pi/2, 20), linspace(0, 0.9, 5));
%! u = phi(:).'; m = sin(alpha(:).').^2; c = cv(:).';
%! Pi_s = elliptic3(u, m, c);
%! elliptic_config('gpu', true);
%! Pi_g = elliptic3(u, m, c);
%! elliptic_config('gpu', false);
%! assert(max(abs(Pi_g - Pi_s)) < 1e-12, 'elliptic3: GPU/serial mismatch');

%!test
%! % ellipj: GPU output must match serial CPU output
%! clear
%! elliptic_config('gpu', false);
%! if ~gpu_available_for_test(), disp('SKIP: no GPU'); return; end
%! [phi, alpha] = meshgrid(linspace(0, 10, 40), linspace(0, pi/2, 40));
%! u = phi(:).'; m = sin(alpha(:).').^2;
%! [Sn_s, Cn_s, Dn_s, Am_s] = ellipj(u, m);
%! elliptic_config('gpu', true);
%! [Sn_g, Cn_g, Dn_g, Am_g] = ellipj(u, m);
%! elliptic_config('gpu', false);
%! assert(max(abs(Sn_g - Sn_s)) < 1e-12, 'ellipj Sn: GPU/serial mismatch');
%! assert(max(abs(Cn_g - Cn_s)) < 1e-12, 'ellipj Cn: GPU/serial mismatch');
%! assert(max(abs(Dn_g - Dn_s)) < 1e-12, 'ellipj Dn: GPU/serial mismatch');
%! assert(max(abs(Am_g - Am_s)) < 1e-12, 'ellipj Am: GPU/serial mismatch');

%!test
%! % jacobiThetaEta: GPU output must match serial CPU output
%! clear
%! elliptic_config('gpu', false);
%! if ~gpu_available_for_test(), disp('SKIP: no GPU'); return; end
%! [phi, qv] = meshgrid(linspace(0.1, pi, 40), linspace(0.01, 0.9, 40));
%! u = phi(:).'; m = qv(:).'.^2;
%! [Th_s, H_s] = jacobiThetaEta(u, m);
%! elliptic_config('gpu', true);
%! [Th_g, H_g] = jacobiThetaEta(u, m);
%! elliptic_config('gpu', false);
%! assert(max(abs(Th_g - Th_s)) < 1e-12, 'jacobiThetaEta Th: GPU/serial mismatch');
%! assert(max(abs(H_g  - H_s )) < 1e-12, 'jacobiThetaEta H:  GPU/serial mismatch');

%!test
%! % Verify GPU is disabled by default and has_gpu() returns false
%! clear
%! elliptic_config('gpu', false);
%! assert(~has_gpu(), 'has_gpu() must return false when gpu config is disabled');
%! assert(~elliptic_config('gpu'), 'gpu config must default to false after clear');
