% BENCH_GPU  Full three-way benchmark: serial CPU / parallel CPU / GPU
%
% Measures wall-clock time, throughput, and hardware utilisation for
% every elliptic function at a wide range of input sizes (10k – 4M pts).
%
% Usage (Octave):
%   pkg load parallel ocl    % load once per session
%   addpath src
%   bench_gpu                % runs everything, writes bench_gpu_results.csv
%
% Usage (MATLAB with Parallel Computing Toolbox + CUDA GPU):
%   parpool(maxNumCompThreads)
%   addpath src
%   bench_gpu
%
% Output:
%   bench_gpu_results.csv   – raw timing + utilisation data

function bench_gpu()
    addpath(fullfile(pwd, 'src'));

    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('  Elliptic Functions — Serial / Parallel / GPU Bench  \n');
    fprintf('=======================================================\n');
    print_system_info();
    fprintf('\n');

    % ── preserve original config ─────────────────────────────────────────
    orig_par   = elliptic_config('parallel');
    orig_gpu   = elliptic_config('gpu');
    orig_chunk = elliptic_config('chunk_size');
    cleanup = onCleanup(@() restore_config(orig_par, orig_gpu, orig_chunk));

    % ── probe hardware ───────────────────────────────────────────────────
    elliptic_config('parallel', true);
    nw = get_nworkers();
    elliptic_config('parallel', false);

    elliptic_config('gpu', true);
    gpu_avail = has_gpu();
    elliptic_config('gpu', false);

    fprintf('Parallel workers : %d\n', nw);
    fprintf('GPU available    : %s\n\n', bool2str(gpu_avail));

    % ── sizes to test (points) ───────────────────────────────────────────
    N_list = [1e4, 5e4, 1e5, 5e5, 1e6, 2e6, 4e6];

    % ── run each function ────────────────────────────────────────────────
    results = {};
    results = [results, bench_fn('elliptic12',     @gen_elliptic12,     N_list, nw, gpu_avail)];
    results = [results, bench_fn('ellipj',         @gen_ellipj,         N_list, nw, gpu_avail)];
    results = [results, bench_fn('elliptic3',      @gen_elliptic3,      N_list, nw, gpu_avail)];
    results = [results, bench_fn('jacobiThetaEta', @gen_jacobiThetaEta, N_list, nw, gpu_avail)];

    print_summary(results, nw, gpu_avail);
    write_csv(results, 'bench_gpu_results.csv');
    fprintf('\nResults saved → bench_gpu_results.csv\n');
end

% ─────────────────────────────────────────────────────────────────────────
%  System info block
% ─────────────────────────────────────────────────────────────────────────
function print_system_info()
    if exist('OCTAVE_VERSION', 'builtin')
        fprintf('Runtime : Octave %s\n', version());
    else
        fprintf('Runtime : MATLAB %s\n', version());
    end
    fprintf('Date    : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    if isunix()
        [~, cpu_str] = system("grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs");
        [~, cores]   = system('nproc');
        fprintf('CPU     : %s (%s cores)\n', strtrim(cpu_str), strtrim(cores));
        [~, mem_str] = system("free -h | awk 'NR==2{print $2}'");
        fprintf('RAM     : %s total\n', strtrim(mem_str));
    end

    [ok, gpu_str] = system('nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv,noheader 2>/dev/null');
    if ok == 0
        parts = strsplit(strtrim(gpu_str), ', ');
        if numel(parts) >= 4
            fprintf('GPU     : %s  VRAM %s  driver %s  compute %s\n', ...
                strtrim(parts{1}), strtrim(parts{2}), strtrim(parts{3}), strtrim(parts{4}));
        end
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  Input generators
% ─────────────────────────────────────────────────────────────────────────
function args = gen_elliptic12(N)
    args = {linspace(0.001, pi/2-0.001, N), linspace(0.001, 0.999, N)};
end
function args = gen_ellipj(N)
    args = {linspace(0.001, 10, N), linspace(0.001, 0.999, N)};
end
function args = gen_elliptic3(N)
    args = {linspace(0.001, pi/2-0.001, N), linspace(0.001, 0.999, N), linspace(0.001, 0.899, N)};
end
function args = gen_jacobiThetaEta(N)
    args = {linspace(0.001, 2, N), linspace(0.001, 0.999, N)};
end

% ─────────────────────────────────────────────────────────────────────────
%  Core benchmark runner for one function
% ─────────────────────────────────────────────────────────────────────────
function rows = bench_fn(fname, gen_fn, N_list, nw, gpu_avail)
    NREPS  = 3;
    WARMUP = 1;

    fprintf('\n─── %s ───\n', fname);
    hdr = '  %-10s  %9s';
    if nw > 0,      hdr = [hdr '  %9s  %7s']; end
    if gpu_avail,   hdr = [hdr '  %9s  %7s  %9s']; end
    fprintf([hdr '\n'], 'N', 'serial(s)', 'par(s)', 'par/ser', 'gpu(s)', 'gpu/ser', 'Mpts/s');
    fprintf('  %s\n', repmat('-', 1, 10 + 11 + 11*nw + 28*gpu_avail));

    rows = {};

    % measure utilisation once at a meaningful size (1M or largest available)
    util_N = N_list(find(N_list >= 1e6, 1));
    if isempty(util_N), util_N = N_list(end); end
    util_done = false;
    util_r = struct('gpu_util', 0, 'gpu_mem_used', 0, 'gpu_mem_total', 0, ...
                    'cpu_ser_pct', 0, 'cpu_par_pct', 0, 'cpu_gpu_pct', 0, ...
                    'par_mem_mb', 0, 'gpu_mem_delta_mib', 0);

    for ni = 1:length(N_list)
        N   = N_list(ni);
        args = gen_fn(N);
        fn   = str2func(fname);

        % serial
        elliptic_config('parallel', false);
        elliptic_config('gpu', false);
        t_serial = timeit_min(fn, args, WARMUP, NREPS);

        % parallel
        t_par = NaN;
        if nw > 0
            chunk = max(floor(N / nw), 1000);
            elliptic_config('parallel', true);
            elliptic_config('gpu', false);
            elliptic_config('chunk_size', chunk);
            t_par = timeit_min(fn, args, WARMUP, NREPS);
            elliptic_config('parallel', false);
        end

        % gpu
        t_gpu = NaN;
        if gpu_avail
            elliptic_config('gpu', true);
            elliptic_config('parallel', false);
            t_gpu = timeit_min(fn, args, WARMUP, NREPS);
            elliptic_config('gpu', false);
        end

        % utilisation snapshot (once at util_N)
        if gpu_avail && ~util_done && N == util_N
            util_r = measure_utilisation(fn, gen_fn, util_N, nw, fname);
            util_done = true;
        end

        sp_par = t_serial / max(t_par, 1e-9);
        sp_gpu = t_serial / max(t_gpu, 1e-9);
        mpts   = N / max(t_gpu, 1e-9) / 1e6;

        line = sprintf('  N=%-9.0f  %9.4f', N, t_serial);
        if nw > 0
            line = [line sprintf('  %9.4f  %6.1fx', t_par, sp_par)];
        end
        if gpu_avail
            line = [line sprintf('  %9.4f  %6.1fx  %9.1f', t_gpu, sp_gpu, mpts)];
        end
        fprintf('%s\n', line);

        r.fname       = fname;
        r.N           = N;
        r.t_serial    = t_serial;
        r.t_par       = t_par;
        r.t_gpu       = t_gpu;
        r.speedup_par = sp_par;
        r.speedup_gpu = sp_gpu;
        r.mpts_gpu    = mpts;
        rows{end+1} = r;
    end

    % print utilisation block
    if gpu_avail
        fprintf('\n  Utilisation @ N=%.0f:\n', util_N);
        fprintf('    Serial  CPU: %3.0f%%  RAM delta: %+.0f MB\n', ...
            util_r.cpu_ser_pct, util_r.par_mem_mb);
        if nw > 0
            fprintf('    Parallel CPU: %3.0f%%  (%d workers)\n', util_r.cpu_par_pct, nw);
        end
        fprintf('    GPU compute: %3d%%  GPU mem used: %d MiB / %d MiB\n', ...
            util_r.gpu_util, util_r.gpu_mem_used, util_r.gpu_mem_total);
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  Timing
% ─────────────────────────────────────────────────────────────────────────
function t = timeit_min(fn, args, warmup, nreps)
    for k = 1:warmup, fn(args{:}); end
    times = zeros(1, nreps);
    for k = 1:nreps
        t0 = tic; fn(args{:}); times(k) = toc(t0);
    end
    t = min(times);
end

% ─────────────────────────────────────────────────────────────────────────
%  Hardware utilisation measurement at a specific N
%  Polls nvidia-smi + /proc/stat while computation runs
% ─────────────────────────────────────────────────────────────────────────
function r = measure_utilisation(fn, gen_fn, N, nw, fname)
    r = struct('gpu_util', 0, 'gpu_mem_used', 0, 'gpu_mem_total', 0, ...
               'cpu_ser_pct', 0, 'cpu_par_pct', 0, 'cpu_gpu_pct', 0, ...
               'par_mem_mb', 0, 'gpu_mem_delta_mib', 0);
    args = gen_fn(N);

    % ── Serial CPU% ──────────────────────────────────────────────────────
    elliptic_config('parallel', false); elliptic_config('gpu', false);
    ct0 = cpu_ticks();
    fn(args{:});
    ct1 = cpu_ticks();
    r.cpu_ser_pct = ticks_to_pct(ct0, ct1);

    % ── Parallel CPU% ────────────────────────────────────────────────────
    if nw > 0
        chunk = max(floor(N/nw), 1000);
        elliptic_config('parallel', true);
        elliptic_config('chunk_size', chunk);
        elliptic_config('gpu', false);
        ct0 = cpu_ticks();
        fn(args{:});
        ct1 = cpu_ticks();
        r.cpu_par_pct = ticks_to_pct(ct0, ct1);
        elliptic_config('parallel', false);
    end

    % ── GPU util + mem ────────────────────────────────────────────────────
    elliptic_config('gpu', true); elliptic_config('parallel', false);

    % Start background nvidia-smi poller (100 ms interval)
    pfile = [tempname(), '.csv'];
    poll_cmd = sprintf( ...
        'bash -c ''while true; do nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits >> "%s" 2>/dev/null; sleep 0.1; done'' &', ...
        pfile);
    system(poll_cmd);

    % Warm up GPU memory (one call before measurement)
    fn(args{:});
    mem_before = gpu_mem_used_mib();

    ct0 = cpu_ticks();
    % Run multiple times so poller gets samples
    for rep = 1:5
        fn(args{:});
    end
    ct1 = cpu_ticks();

    mem_after = gpu_mem_used_mib();
    r.gpu_mem_delta_mib = mem_after - mem_before;
    r.cpu_gpu_pct = ticks_to_pct(ct0, ct1);

    % Kill poller
    system('pkill -f "nvidia-smi --query-gpu=utilization.gpu" 2>/dev/null; true');
    pause(0.2);

    elliptic_config('gpu', false);

    % Parse poll file
    if exist(pfile, 'file')
        fid = fopen(pfile, 'r');
        raw = textscan(fid, '%f,%f,%f', 'CollectOutput', true);
        fclose(fid); delete(pfile);
        if ~isempty(raw{1}) && size(raw{1}, 1) > 0
            data = raw{1};
            % Exclude idle (first/last) samples, use 80th percentile for util
            if size(data,1) > 4
                data = data(2:end-2, :);
            end
            r.gpu_util      = round(prctile(data(:,1), 80));
            r.gpu_mem_used  = round(max(data(:,2)));
            r.gpu_mem_total = round(data(1,3));
        end
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  CPU tick helpers
% ─────────────────────────────────────────────────────────────────────────
function t = cpu_ticks()
    t = [];
    if ~isunix(), return; end
    fid = fopen('/proc/stat', 'r');
    if fid < 0, return; end
    line = fgetl(fid); fclose(fid);
    parts = strsplit(strtrim(line));
    if numel(parts) >= 5
        t = cellfun(@str2double, parts(2:end));
    end
end

function pct = ticks_to_pct(t0, t1)
    pct = 0;
    if isempty(t0) || isempty(t1) || numel(t0) ~= numel(t1), return; end
    delta = t1 - t0;
    total = sum(delta);
    if total <= 0, return; end
    pct = 100 * (1 - delta(4) / total);   % col 4 = idle
end

% ─────────────────────────────────────────────────────────────────────────
%  GPU memory query
% ─────────────────────────────────────────────────────────────────────────
function mib = gpu_mem_used_mib()
    mib = 0;
    [ok, s] = system('nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null');
    if ok == 0
        mib = str2double(strtrim(s));
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  Summary table
% ─────────────────────────────────────────────────────────────────────────
function print_summary(results, nw, gpu_avail)
    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('  SUMMARY — Best speedups at N >= 500k\n');
    fprintf('=======================================================\n');

    hdr = '  %-20s  %10s';
    if nw > 0,    hdr = [hdr '  %10s']; end
    if gpu_avail, hdr = [hdr '  %10s  %12s']; end
    args = {'Function', 'serial(s)'};
    if nw > 0,    args = [args, {'par speedup'}]; end
    if gpu_avail, args = [args, {'gpu speedup', 'GPU Mpts/s'}]; end
    fprintf([hdr '\n'], args{:});
    fprintf('  %s\n', repmat('-', 1, 22 + 12 + 12*nw + 24*gpu_avail));

    fnames = unique(cellfun(@(r) r.fname, results, 'UniformOutput', false));
    for fi = 1:numel(fnames)
        fn   = fnames{fi};
        mask = cellfun(@(r) strcmp(r.fname, fn) && r.N >= 5e5, results);
        sub  = results(mask);
        if isempty(sub), continue; end

        t_ser  = min(cellfun(@(r) r.t_serial,    sub));
        sp_par = max(cellfun(@(r) r.speedup_par, sub));
        sp_gpu = max(cellfun(@(r) r.speedup_gpu, sub));
        mpts   = max(cellfun(@(r) r.mpts_gpu,    sub));

        line = sprintf('  %-20s  %9.3fs', fn, t_ser);
        if nw > 0
            line = [line sprintf('  %9.1fx', sp_par)];
        end
        if gpu_avail
            if isnan(sp_gpu)
                line = [line sprintf('  %10s  %12s', 'N/A', 'N/A')];
            else
                line = [line sprintf('  %9.1fx  %10.1f M', sp_gpu, mpts)];
            end
        end
        fprintf('%s\n', line);
    end
    fprintf('\n');
    fprintf('  serial(s)  = fastest serial run at N >= 500k\n');
    fprintf('  par speedup = serial / parallel  (best-of-3 min wall-clock)\n');
    if gpu_avail
        fprintf('  gpu speedup = serial / GPU       (includes host↔device transfer)\n');
        fprintf('  GPU Mpts/s  = million points per second through GPU path\n');
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  CSV export
% ─────────────────────────────────────────────────────────────────────────
function write_csv(results, outfile)
    fid = fopen(outfile, 'w');
    fprintf(fid, ['function,N,t_serial_s,t_parallel_s,t_gpu_s,' ...
        'speedup_par,speedup_gpu,mpts_per_s_gpu\n']);
    for k = 1:numel(results)
        r = results{k};
        fprintf(fid, '%s,%d,%.6f,%.6f,%.6f,%.4f,%.4f,%.2f\n', ...
            r.fname, r.N, r.t_serial, r.t_par, r.t_gpu, ...
            r.speedup_par, r.speedup_gpu, r.mpts_gpu);
    end
    fclose(fid);
end

% ─────────────────────────────────────────────────────────────────────────
%  Utilities
% ─────────────────────────────────────────────────────────────────────────
function s = bool2str(b)
    if b, s = 'yes'; else, s = 'no'; end
end

function p = prctile(v, p_pct)
    % Simple percentile (avoids Statistics Toolbox dependency)
    sv = sort(v(:));
    idx = max(1, round(p_pct/100 * numel(sv)));
    p = sv(idx);
end

function restore_config(par, gpu, chunk)
    elliptic_config('parallel', par);
    elliptic_config('gpu', gpu);
    elliptic_config('chunk_size', chunk);
end
