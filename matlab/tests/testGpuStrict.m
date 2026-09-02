function testGpuStrict()
%TESTGPUSTRICT  Every GPU code path under a strict device-array stub.
%   See gpu_stub/README.md.  Each function is evaluated on the CPU path and
%   on the GPU path (with the stub) on the same inputs, including the
%   large-u / m -> 1 / tiny-m cases that broke earlier GPU kernels, and the
%   two must agree to 1e-13 relative (they are bit-identical for most).
end

%!test
%! % locate the stub from the source directory: mfilename is empty inside test
%! % blocks when test() is called with a full path (CI), so relative paths fail
%! src = fileparts(which('elliptic12'));
%! here = fullfile(src, '..', 'tests');
%! addpath(fullfile(here, 'gpu_stub'), '-begin');
%! unwind_protect
%!   rand('seed', 5); N = 300;
%!   u = rand(1,N)*20 - 10; m = rand(1,N)*(1-2e-6) + 1e-6; n = rand(1,N)*0.9; z = rand(1,N)*3;
%!   u = [u, 1e6+0.123, 9.375, pi/2-1e-9, -4.9, 0, 1000000.123]; m = [m, 1-eps/2, 1-eps/2, 1-eps/2, 1e-12, 0.5, 3e-16];
%!   n = [n 0.3 0.3 0.3 0.3 0.3 0.3]; z = [z 0.7 0.7 0.7 0.7 0.7 0.7];
%!   tests = {
%!     'elliptic12',  @() nthargout(1:3, @elliptic12, u, m);
%!     'ellipj',      @() nthargout(1:4, @ellipj, u, m);
%!     'elliptic3',   @() elliptic3(u, m, n);
%!     'ellipticBDJ', @() nthargout(1:3, @ellipticBDJ, u, m, n);
%!     'ellipticBD',  @() nthargout(1:3, @ellipticBD, m);
%!     'jacobiEDJ',   @() nthargout(1:3, @jacobiEDJ, u, m, n);
%!     'theta',       @() theta(1, u, m);
%!     'theta_prime', @() nthargout(1:2, @theta_prime, 2, u, m);
%!     'jacobiThetaEta', @() nthargout(1:2, @jacobiThetaEta, u, m);
%!     'nomeq',       @() nomeq(m);
%!     'inversenomeq',@() inversenomeq(m*0.7);
%!     'elliptic12i', @() nthargout(1:3, @elliptic12i, u + 1i*z, m);
%!     'ellipji',     @() nthargout(1:3, @ellipji, u + 1i*z, m);
%!     'weierstrassP',@() weierstrassP(z, 1.5, -0.25, -1.25);
%!     'weierstrassZeta', @() weierstrassZeta(z, 1.5, -0.25, -1.25);
%!     'weierstrassSigma', @() weierstrassSigma(z, 1.5, -0.25, -1.25);
%!     'weierstrassPPrime', @() weierstrassPPrime(z, 1.5, -0.25, -1.25);
%!     'inverselliptic2', @() inverselliptic2(z, m);
%!     'cel',         @() cel(sqrt(1-m), n, 1, 0.5);
%!     'arclength_ellipse', @() arclength_ellipse(z+0.1, z+0.5, u, z);
%!     'elliptic12 NaN', @() nthargout(1:3, @elliptic12, [0.3 0.5 0.7 NaN], [0.2 NaN 0.4 0.5]);
%!     'ellipj NaN',     @() nthargout(1:4, @ellipj, [0.3 0.5 0.7 NaN], [0.2 NaN 0.4 0.5]);
%!     'theta NaN',      @() theta(1, [0.3 0.5], [0.2 NaN]);
%!     'elliptic3 NaN',  @() elliptic3([0.3 0.5], [0.2 NaN], 0.3);
%!     'elliptic3 c<0 / near pole', @() elliptic3([1 1 1 4 1.5707 1.2 0.4], [0.5 0.5 0.5 0.9 1-1e-9 0.999999 0.3], [-0.5 -3 -100 -100 0.3 0.999999 0.5]);
%!   };
%!   for t = 1:rows(tests)
%!     name = tests{t,1}; f = tests{t,2};
%!     elliptic_config('gpu', false); ref = f();
%!     elliptic_config('gpu', true);  got = f();
%!     elliptic_config('gpu', false);
%!     if ~iscell(ref), ref = {ref}; got = {got}; end
%!     for k = 1:numel(ref)
%!       g = got{k};
%!       assert(~isa(g, 'gpuArray'), sprintf('%s: output %d is still a device array', name, k));
%!       r = ref{k};
%!       assert(isequal(isfinite(r), isfinite(g)), sprintf('%s: finite pattern differs on the GPU path', name));
%!       ok = isfinite(r);
%!       d = max([0; abs(g(ok)(:) - r(ok)(:)) ./ max(1, abs(r(ok)(:)))]);
%!       assert(d < 1e-13, sprintf('%s: GPU path differs from CPU by %.2e', name, d));
%!     end
%!   end
%! unwind_protect_cleanup
%!   elliptic_config('gpu', false);
%!   rmpath(fullfile(here, 'gpu_stub'));
%! end_unwind_protect

%!test
%! % the stub itself must reject what ocl rejects, or the test above proves nothing
%! here = fullfile(fileparts(which('elliptic12')), '..', 'tests');
%! addpath(fullfile(here, 'gpu_stub'), '-begin');
%! unwind_protect
%!   caught = false;
%!   try, x = gpuArray([1 2 3]) .* [1 2 3]; catch, caught = true; end
%!   assert(caught, 'stub must reject device .* host matrix');
%!   caught = false;
%!   try, g = gpuArray([1 2 3]); y = g(logical([1 0 1])); catch, caught = true; end
%!   assert(caught, 'stub must reject logical indexing');
%!   assert(isequal(gather(gpuArray([1 2 3]) .* 2 + gpuArray([1 1 1])), [3 5 7]));
%! unwind_protect_cleanup
%!   rmpath(fullfile(here, 'gpu_stub'));
%! end_unwind_protect
