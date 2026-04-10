function tf = has_gpu()
%HAS_GPU  Returns true if GPU is enabled in config and a compatible GPU is available.
%
%   TF = HAS_GPU() returns true only when both conditions hold:
%     (a) elliptic_config('gpu') is true  — user must opt in explicitly
%     (b) a compatible GPU is present and accessible
%
%   MATLAB: checks canUseGPU() (R2020b+) or gpuDeviceCount() > 0 as fallback.
%   Octave: attempts to load the 'ocl' package (OpenCL) and probe the GPU.
%           Falls back to CPU silently if ocl is unavailable.
%
%   Enable GPU mode:
%       elliptic_config('gpu', true);
%
%   See also ELLIPTIC_CONFIG, GET_NWORKERS.

    tf = false;

    if ~elliptic_config('gpu')
        return;
    end

    if exist('OCTAVE_VERSION', 'builtin')
        % Octave: try the 'ocl' Forge package (OpenCL via gpuArray/gather API)
        try
            pkg('load', 'ocl');
            tmp = gpuArray(1.0);
            gather(tmp);   % confirm round-trip works
            tf = true;
        catch
            % ocl not installed or no OpenCL device — fall back to CPU silently
        end
        return;
    end

    % MATLAB: prefer canUseGPU (R2020b+), fall back to gpuDeviceCount
    if exist('canUseGPU', 'builtin') || exist('canUseGPU', 'file')
        tf = canUseGPU();
    elseif exist('gpuDeviceCount', 'file')
        tf = gpuDeviceCount() > 0;
    end
