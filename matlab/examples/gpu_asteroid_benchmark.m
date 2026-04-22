function gpu_asteroid_benchmark(varargin)
% GPU_ASTEROID_BENCHMARK  Massive Keplerian propagation + elliptic arc length.
%
%   gpu_asteroid_benchmark()
%   gpu_asteroid_benchmark('bodies', [1e4 1e5 1e6], 'epochs', [128 512 2048])
%
% Loads an asteroid orbital-element dataset (npz produced by
% scripts/prepare_asteroid_dataset.py), propagates each body on a time
% grid via Kepler's equation, and computes orbital arc lengths using
% elliptic12. Runs both on the CPU and on a CUDA gpuArray, saving
% timings to examples/gpu-asteroid-swarm/data/benchmark_results_matlab.json
%
% This is an HONEST throughput demo: it is a batched pure-Keplerian
% propagator, not a full n-body integrator.

  p = inputParser;
  addParameter(p, 'bodies', [1e4 1e5 1e6]);
  addParameter(p, 'epochs', [128 512 2048]);
  addParameter(p, 'horizonYears', 10);
  addParameter(p, 'npz', fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
               'examples', 'gpu-asteroid-swarm', 'data', 'asteroids_full.npz'));
  addParameter(p, 'out', fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
               'examples', 'gpu-asteroid-swarm', 'data', 'benchmark_results_matlab.json'));
  parse(p, varargin{:});
  opts = p.Results;

  addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

  hasGPU = exist('gpuDeviceCount','file') && gpuDeviceCount > 0;
  if hasGPU
    dev = gpuDevice();
    fprintf('GPU: %s (%.1f GB)\n', dev.Name, dev.TotalMemory / 2^30);
  else
    fprintf('No GPU detected — running CPU-only benchmark\n');
  end

  [a_all, e_all, i_all, Om_all, w_all, M0_all, epoch_all] = load_elements(opts.npz);
  fprintf('Loaded %d bodies\n', numel(a_all));

  tTemplate = linspace(0, opts.horizonYears*365.25, max(opts.epochs));
  rows = struct([]);

  for N = opts.bodies(:).'
    N = min(N, numel(a_all));
    idx = randperm(numel(a_all), N);
    a  = a_all(idx);  e  = e_all(idx);  inc = i_all(idx);
    Om = Om_all(idx); w  = w_all(idx);  M0 = M0_all(idx);
    ep = epoch_all(idx);

    for T = opts.epochs(:).'
      t = tTemplate(1:T);

      % CPU
      tic;
      [xc, yc, zc, Ec] = propagate(a, e, inc, Om, w, M0, ep, t);
      Lc = arclen(a, e, Ec);
      t_cpu = toc;

      % GPU
      t_gpu = NaN;
      if hasGPU
        gd = @(v) gpuArray(double(v));
        wait(dev);
        tic;
        [xg, yg, zg, Eg] = propagate(gd(a), gd(e), gd(inc), gd(Om), ...
                                     gd(w), gd(M0), gd(ep), gd(t));
        Lg = arclen(gd(a), gd(e), Eg);
        wait(dev);
        t_gpu = toc;
      end

      row = struct( ...
        'bodies', N, 'epochs', T, 'body_epochs', N*T, ...
        't_cpu_ms', 1e3*t_cpu, 't_gpu_ms', 1e3*t_gpu, ...
        'throughput_meps_cpu', (N*T)/t_cpu/1e6, ...
        'throughput_meps_gpu', (N*T)/t_gpu/1e6, ...
        'speedup', t_cpu/t_gpu);
      rows = [rows, row];
      fprintf('  N=%-8d T=%-5d  CPU=%8.1f ms  GPU=%8.1f ms  speedup=%5.1fx\n', ...
              N, T, 1e3*t_cpu, 1e3*t_gpu, t_cpu/t_gpu);
    end
  end

  payload = struct( ...
    'matlab_version', version, ...
    'has_gpu', hasGPU, ...
    'gpu_name', ternary(hasGPU, dev.Name, ''), ...
    'results', rows);
  fid = fopen(opts.out, 'w');
  fprintf(fid, '%s', jsonencode(payload, 'PrettyPrint', true));
  fclose(fid);
  fprintf('Wrote %s\n', opts.out);
end

% -------------------------------------------------------------------------
function [x, y, z, E] = propagate(a, e, inc, Om, w, M0, epoch, t)
  GM_SUN = 0.0002959122082855911;  % AU^3/day^2
  a = a(:); e = e(:); inc = inc(:); Om = Om(:); w = w(:); M0 = M0(:); epoch = epoch(:);
  t = t(:).';                      % row

  n = sqrt(GM_SUN ./ a.^3);        % (N,1)
  M = M0 + n .* (t - epoch);       % (N,T)
  M = mod(M + pi, 2*pi) - pi;

  E = M + e .* sin(M);
  for k = 1:12
    E = E - (E - e.*sin(E) - M) ./ (1 - e.*cos(E));
  end
  ca = cos(E); sa = sin(E);
  b  = a .* sqrt(1 - e.^2);
  x_op = a .* (ca - e);
  y_op = b .* sa;

  cw = cos(w); sw = sin(w);
  cO = cos(Om); sO = sin(Om);
  ci = cos(inc); si = sin(inc);
  R11 =  cO.*cw - sO.*sw.*ci;   R12 = -cO.*sw - sO.*cw.*ci;
  R21 =  sO.*cw + cO.*sw.*ci;   R22 = -sO.*sw + cO.*cw.*ci;
  R31 =  sw.*si;                R32 =  cw.*si;
  x = R11.*x_op + R12.*y_op;
  y = R21.*x_op + R22.*y_op;
  z = R31.*x_op + R32.*y_op;
end

function L = arclen(a, e, E)
  m = (e.^2) .* ones(size(E));
  [~, Ephi] = elliptic12(E, m);
  L = a .* Ephi;
end

% -------------------------------------------------------------------------
function [a, e, inc, Om, w, M0, epoch] = load_elements(npz)
  if exist(npz, 'file')
    % MATLAB can't read npz natively; fall back to per-platform npz reader
    % For robustness we allow a CSV sibling file produced by the Python prep.
    csvPath = strrep(npz, '.npz', '.csv');
    if exist(csvPath, 'file')
      T = readtable(csvPath);
      a = T.a; e = T.e; inc = T.i; Om = T.om; w = T.w; M0 = T.ma; epoch = T.epoch;
      return;
    end
  end
  % synthetic fallback
  rng(20260421);
  N = 1e6;
  a  = 1.5 + 3.0 * rand(N, 1);
  e  = min(0.95, betarnd(2, 15, N, 1) * 1.2);
  inc = deg2rad(abs(10 * randn(N, 1)));
  Om  = 2*pi * rand(N, 1);
  w   = 2*pi * rand(N, 1);
  M0  = 2*pi * rand(N, 1);
  epoch = 2460000.5 * ones(N, 1);
end

function v = ternary(c, a, b)
  if c, v = a; else, v = b; end
end
