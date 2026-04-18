% BENCH  Benchmark suite for elliptic functions library
% Runs each function with varying input sizes and reports timing.
% If a parallel pool is active, parallelism kicks in automatically
% for inputs larger than elliptic_config('chunk_size').

function bench()
    addpath(fullfile(pwd, 'src'));
    fprintf('\n=== Elliptic Functions Benchmark ===\n');
    if exist('OCTAVE_VERSION', 'builtin')
        fprintf('Octave %s on %s\n', version, computer);
    else
        fprintf('MATLAB %s on %s\n', version, computer);
    end
    fprintf('Date: %s\n', datestr(now));

    nw = get_nworkers();
    fprintf('Parallel workers: %d\n', nw);
    fprintf('Chunk size: %d\n\n', elliptic_config('chunk_size'));

    bench_elliptic12(nw);
    bench_elliptic3(nw);
    bench_ellipj(nw);
    bench_jacobiThetaEta(nw);

    fprintf('\n=== Benchmark Complete ===\n');
end

function t = timeit_fn(fn, nreps)
    if nargin < 2, nreps = 5; end
    fn();  % warmup
    times = zeros(1, nreps);
    for k = 1:nreps
        tic;
        fn();
        times(k) = toc;
    end
    t = mean(times);
end

function bench_elliptic12(nw)
    fprintf('--- elliptic12 ---\n');
    sizes = [100, 500, 1000, 2000];
    for s = sizes
        [phi, alpha] = meshgrid(linspace(0, pi/2, s), linspace(0, 89, s));
        m = sin(pi/180*alpha).^2;
        u_flat = phi(:).'; m_flat = m(:).';

        t = timeit_fn(@() elliptic12(u_flat, m_flat));
        fprintf('  %5dx%-5d (%7d pts): %.4fs', s, s, s*s, t);
        if nw > 1
            fprintf('  [parallel: %d workers]', nw);
        end
        fprintf('\n');
    end
end

function bench_elliptic3(nw)
    fprintf('\n--- elliptic3 ---\n');
    sizes = [100, 1000, 10000, 100000];
    for s = sizes
        u = linspace(0.01, pi/2 - 0.01, s);
        m = linspace(0.01, 0.99, s);
        c = linspace(0.01, 0.99, s);

        t = timeit_fn(@() elliptic3(u, m, c));
        fprintf('  %7d pts: %.6fs', s, t);
        if nw > 1
            fprintf('  [parallel: %d workers]', nw);
        end
        fprintf('\n');
    end
end

function bench_ellipj(nw)
    fprintf('\n--- ellipj ---\n');
    sizes = [100, 500, 1000, 2000];
    for s = sizes
        [phi, alpha] = meshgrid(linspace(0, 4, s), linspace(0, 89, s));
        m = sin(pi/180*alpha).^2;
        u = phi(:).'; m = m(:).';

        t = timeit_fn(@() ellipj(u, m));
        fprintf('  %5dx%-5d (%7d pts): %.4fs', s, s, s*s, t);
        if nw > 1
            fprintf('  [parallel: %d workers]', nw);
        end
        fprintf('\n');
    end
end

function bench_jacobiThetaEta(nw)
    fprintf('\n--- jacobiThetaEta ---\n');
    sizes = [100, 500, 1000, 2000];
    for s = sizes
        [uu, mm] = meshgrid(linspace(0.01, 1, s), linspace(0.01, 0.99, s));
        u = uu(:).'; m = mm(:).';

        t = timeit_fn(@() jacobiThetaEta(u, m));
        fprintf('  %5dx%-5d (%7d pts): %.4fs', s, s, s*s, t);
        if nw > 1
            fprintf('  [parallel: %d workers]', nw);
        end
        fprintf('\n');
    end
end
