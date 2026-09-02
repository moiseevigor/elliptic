function r = has_gpu()
%HAS_GPU  Test stub: the strict device stub is "available" whenever the
%   configuration flag is on (the real has_gpu also needs the ocl package).
r = elliptic_config('gpu');
end
