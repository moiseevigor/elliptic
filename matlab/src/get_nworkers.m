function n = get_nworkers()
%GET_NWORKERS  Returns the number of available parallel workers.
%   N = GET_NWORKERS() returns the number of workers in the current
%   parallel pool (MATLAB) or the number of available CPU cores (Octave).
%   Returns 0 if no parallel support is available or not enabled.
%
%   Parallelism must be enabled via elliptic_config('parallel', true).
%   Without this, always returns 0 (serial mode).
%
%   MATLAB: Requires Parallel Computing Toolbox with an active pool.
%           Use `parpool(N)` to start a pool before calling parallel functions.
%
%   Octave: Requires the `parallel` package (`pkg install -forge parallel`).
%           Returns nproc() (number of available cores).
%
%   See also ELLIPTIC_CONFIG, HAS_GPU, PARPOOL.

    n = 0;

    if ~elliptic_config('parallel')
        return;
    end

    if exist('OCTAVE_VERSION', 'builtin')
        % Octave: check for parallel package
        try
            pkg('load', 'parallel');
            n = nproc();
        catch
            n = 0;
        end
    else
        % MATLAB: check for Parallel Computing Toolbox
        if exist('gcp', 'file')
            pool = gcp('nocreate');
            if ~isempty(pool)
                n = pool.NumWorkers;
            end
        end
    end
