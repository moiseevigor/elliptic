% BENCH  Benchmark suite for elliptic functions library
% Runs each function with varying input sizes and reports timing.

function bench()
    addpath(fullfile(pwd, 'src'));
    fprintf('\n=== Elliptic Functions Benchmark ===\n');
    fprintf('Octave %s on %s\n', version, computer);
    fprintf('Date: %s\n\n', datestr(now));

    % --- elliptic12 ---
    fprintf('--- elliptic12 ---\n');
    sizes = [100, 500, 1000, 2000];
    for s = sizes
        [phi, alpha] = meshgrid(linspace(0, pi/2, s), linspace(0, 89, s));
        m = sin(pi/180*alpha).^2;
        % warmup
        [F,E,Z] = elliptic12(phi(:).', m(:).');
        times = zeros(1, 5);
        for k = 1:5
            tic;
            [F,E,Z] = elliptic12(phi(:).', m(:).');
            times(k) = toc;
        end
        fprintf('  %5dx%-5d (%7d pts): mean=%.4fs  std=%.4fs\n', s, s, s*s, mean(times), std(times));
    end

    % --- elliptic3 ---
    fprintf('\n--- elliptic3 ---\n');
    sizes = [100, 1000, 10000, 100000];
    for s = sizes
        u = linspace(0.01, pi/2 - 0.01, s);
        m = linspace(0.01, 0.99, s);
        c = linspace(0.01, 0.99, s);
        % warmup
        Pi = elliptic3(u, m, c);
        times = zeros(1, 5);
        for k = 1:5
            tic;
            Pi = elliptic3(u, m, c);
            times(k) = toc;
        end
        fprintf('  %5d pts: mean=%.6fs  std=%.6fs\n', s, mean(times), std(times));
    end

    % --- ellipj ---
    fprintf('\n--- ellipj ---\n');
    sizes = [100, 500, 1000, 2000];
    for s = sizes
        [phi, alpha] = meshgrid(linspace(0, 4, s), linspace(0, 89, s));
        m = sin(pi/180*alpha).^2;
        u = phi(:).';
        m = m(:).';
        % warmup
        [Sn,Cn,Dn,Am] = ellipj(u, m);
        times = zeros(1, 5);
        for k = 1:5
            tic;
            [Sn,Cn,Dn,Am] = ellipj(u, m);
            times(k) = toc;
        end
        fprintf('  %5dx%-5d (%7d pts): mean=%.4fs  std=%.4fs\n', s, s, s*s, mean(times), std(times));
    end

    % --- jacobiThetaEta ---
    fprintf('\n--- jacobiThetaEta ---\n');
    sizes = [100, 500, 1000, 2000];
    for s = sizes
        [uu, mm] = meshgrid(linspace(0.01, 1, s), linspace(0.01, 0.99, s));
        u = uu(:).';
        m = mm(:).';
        % warmup
        [Th, H] = jacobiThetaEta(u, m);
        times = zeros(1, 5);
        for k = 1:5
            tic;
            [Th, H] = jacobiThetaEta(u, m);
            times(k) = toc;
        end
        fprintf('  %5dx%-5d (%7d pts): mean=%.4fs  std=%.4fs\n', s, s, s*s, mean(times), std(times));
    end

    fprintf('\n=== Benchmark Complete ===\n');
end
