function val = elliptic_config(key, value)
%ELLIPTIC_CONFIG  Get/set library configuration for parallel and GPU modes.
%   cfg = ELLIPTIC_CONFIG() returns the current configuration struct.
%
%   val = ELLIPTIC_CONFIG(KEY) returns the value of a configuration key.
%
%   ELLIPTIC_CONFIG(KEY, VALUE) sets a configuration key to a new value.
%
%   Available keys:
%     'parallel'   — logical, enable multi-core via parfor (default: false)
%     'gpu'        — logical, enable GPU acceleration (default: false)
%     'chunk_size' — minimum elements per parallel chunk (default: 10000)
%
%   Example:
%       elliptic_config('parallel', true);   % enable parfor
%       elliptic_config('gpu', true);        % enable gpuArray
%       elliptic_config('parallel')          % query current setting
%
%   See also PAR_ELLIPTIC12, GPU_ELLIPTIC3, GET_NWORKERS, HAS_GPU.

    persistent config;
    if isempty(config)
        config = struct('parallel', false, 'gpu', false, 'chunk_size', 10000);
    end

    if nargin == 0
        val = config;
        return;
    end

    if nargin == 1
        if ~isfield(config, key)
            error('Unknown configuration key: %s', key);
        end
        val = config.(key);
        return;
    end

    if ~isfield(config, key)
        error('Unknown configuration key: %s', key);
    end
    config.(key) = value;
    if nargout > 0
        val = value;
    end
end
