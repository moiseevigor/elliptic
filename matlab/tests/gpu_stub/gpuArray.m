classdef gpuArray
  % Strict stand-in for an ocl/OpenCL device array: elementwise ops work only
  % between gpuArrays or with plain scalars; mixing with a host matrix errors
  % exactly like the ocl package does ("not implemented for 'ocl matrix' by
  % 'matrix'"), logical indexing errors, isreal is false.
  properties
    d
  end
  methods
    function g = gpuArray(x)
      if isa(x, 'gpuArray'), g.d = x.d; else, g.d = double(x); end
    end
    function x = gather(g), x = g.d; end
    function x = double(g), error('gpuArray: double() of a device array; use gather()'); end
    function r = isreal(~), r = false; end
    function r = isa_gpu(~), r = true; end
    function s = size(g, varargin), s = size(g.d, varargin{:}); end
    function n = numel(g), n = numel(g.d); end
    function n = length(g), n = length(g.d); end
    function r = isempty(g), r = isempty(g.d); end
    function r = isscalar(g), r = isscalar(g.d); end
    function r = ndims(g), r = ndims(g.d); end
    function r = columns(g), r = columns(g.d); end
    function r = rows(g), r = rows(g.d); end
    function r = isnumeric(~), r = true; end
    function r = isfloat(~), r = true; end
    function r = islogical(~), r = false; end
    function disp(g), disp('gpuArray (strict stub)'); disp(g.d); end
    function display(g), disp(g); end
    % ---- binary elementwise ----
    function r = plus(a, b),    r = gpuArray(gpuArray.bin(a, b, @plus)); end
    function r = minus(a, b),   r = gpuArray(gpuArray.bin(a, b, @minus)); end
    function r = times(a, b),   r = gpuArray(gpuArray.bin(a, b, @times)); end
    function r = rdivide(a, b), r = gpuArray(gpuArray.bin(a, b, @rdivide)); end
    function r = ldivide(a, b), r = gpuArray(gpuArray.bin(a, b, @ldivide)); end
    function r = power(a, b),   r = gpuArray(gpuArray.bin(a, b, @power)); end
    function r = mtimes(a, b)
      if isscalar(a) || isscalar(b), r = times(a, b); else, error('gpuArray: matrix product not supported by ocl'); end
    end
    function r = mrdivide(a, b)
      if isscalar(b), r = rdivide(a, b); else, error('gpuArray: mrdivide not supported'); end
    end
    function r = mpower(a, b)
      if isscalar(a) && isscalar(b), r = power(a, b); else, error('gpuArray: mpower not supported'); end
    end
    function r = lt(a, b), r = gpuArray(gpuArray.bin(a, b, @lt)); end
    function r = gt(a, b), r = gpuArray(gpuArray.bin(a, b, @gt)); end
    function r = le(a, b), r = gpuArray(gpuArray.bin(a, b, @le)); end
    function r = ge(a, b), r = gpuArray(gpuArray.bin(a, b, @ge)); end
    function r = eq(a, b), r = gpuArray(gpuArray.bin(a, b, @eq)); end
    function r = ne(a, b), r = gpuArray(gpuArray.bin(a, b, @ne)); end
    function r = and(a, b), r = gpuArray(gpuArray.bin(a, b, @and)); end
    function r = or(a, b), r = gpuArray(gpuArray.bin(a, b, @or)); end
    function r = max(a, b, varargin)
      if nargin == 1, r = gpuArray(max(a.d)); else, r = gpuArray(gpuArray.bin(a, b, @max)); end
    end
    function r = min(a, b, varargin)
      if nargin == 1, r = gpuArray(min(a.d)); else, r = gpuArray(gpuArray.bin(a, b, @min)); end
    end
    function r = atan2(a, b), r = gpuArray(gpuArray.bin(a, b, @atan2)); end
    function r = mod(a, b), r = gpuArray(gpuArray.bin(a, b, @mod)); end
    % ---- unary ----
    function r = uminus(g), r = gpuArray(-g.d); end
    function r = uplus(g), r = g; end
    function r = not(g), r = gpuArray(~g.d); end
    function r = sqrt(g), r = gpuArray(sqrt(g.d)); end
    function r = sin(g), r = gpuArray(sin(g.d)); end
    function r = cos(g), r = gpuArray(cos(g.d)); end
    function r = tan(g), r = gpuArray(tan(g.d)); end
    function r = atan(g), r = gpuArray(atan(g.d)); end
    function r = asin(g), r = gpuArray(asin(g.d)); end
    function r = acos(g), r = gpuArray(acos(g.d)); end
    function r = sinh(g), r = gpuArray(sinh(g.d)); end
    function r = cosh(g), r = gpuArray(cosh(g.d)); end
    function r = tanh(g), r = gpuArray(tanh(g.d)); end
    function r = exp(g), r = gpuArray(exp(g.d)); end
    function r = log(g), r = gpuArray(log(g.d)); end
    function r = log1p(g), r = gpuArray(log1p(g.d)); end
    function r = abs(g), r = gpuArray(abs(g.d)); end
    function r = sign(g), r = gpuArray(sign(g.d)); end
    function r = floor(g), r = gpuArray(floor(g.d)); end
    function r = ceil(g), r = gpuArray(ceil(g.d)); end
    function r = round(g), r = gpuArray(round(g.d)); end
    function r = fix(g), r = gpuArray(fix(g.d)); end
    function r = real(g), r = gpuArray(real(g.d)); end
    function r = imag(g), r = gpuArray(imag(g.d)); end
    function r = isnan(g), r = gpuArray(isnan(g.d)); end
    function r = isinf(g), r = gpuArray(isinf(g.d)); end
    function r = isfinite(g), r = gpuArray(isfinite(g.d)); end
    function r = transpose(g), r = gpuArray(g.d.'); end
    function r = ctranspose(g), r = gpuArray(g.d'); end
    function r = sum(g, varargin), r = gpuArray(sum(g.d, varargin{:})); end
    function r = prod(g, varargin), r = gpuArray(prod(g.d, varargin{:})); end
    function r = any(g, varargin), r = any(g.d, varargin{:}); end
    function r = all(g, varargin), r = all(g.d, varargin{:}); end
    function r = reshape(g, varargin), r = gpuArray(reshape(g.d, varargin{:})); end
    function r = repmat(g, varargin), r = gpuArray(repmat(g.d, varargin{:})); end
    function r = horzcat(varargin), r = gpuArray(horzcat(gpuArray.cellu(varargin){:})); end
    function r = vertcat(varargin), r = gpuArray(vertcat(gpuArray.cellu(varargin){:})); end
    function r = cat(dim, varargin), r = gpuArray(cat(dim, gpuArray.cellu(varargin){:})); end
    % ---- indexing: numeric/colon only (ocl has no logical indexing) ----
    function r = subsref(g, s)
      switch s(1).type
        case '()'
          for k = 1:numel(s(1).subs)
            ix = s(1).subs{k};
            if islogical(ix) || isa(ix, 'gpuArray'), error('gpuArray: logical / device-array indexing not supported by ocl'); end
          end
          r = gpuArray(g.d(s(1).subs{:}));
          if numel(s) > 1, r = subsref(r, s(2:end)); end
        case '.'
          if strcmp(s(1).subs, 'd'), r = g.d; else, error('gpuArray: no field %s', s(1).subs); end
        otherwise, error('gpuArray: unsupported indexing');
      end
    end
    function g = subsasgn(g, s, v)
      if ~strcmp(s(1).type, '()'), error('gpuArray: unsupported assignment'); end
      for k = 1:numel(s(1).subs)
        ix = s(1).subs{k};
        if islogical(ix) || isa(ix, 'gpuArray'), error('gpuArray: logical / device-array indexed assignment not supported by ocl'); end
      end
      if isa(v, 'gpuArray'), v = v.d; elseif ~isscalar(v), error('gpuArray: assigning a host matrix into a device array'); end
      g.d(s(1).subs{:}) = v;
    end
    function n = end(g, k, n_), if n_ == 1, n = numel(g.d); else, n = size(g.d, k); end, end
  end
  methods (Static)
    function x = bin(a, b, op)
      ga = isa(a, 'gpuArray'); gb = isa(b, 'gpuArray');
      if ga, av = a.d; else, av = a; end
      if gb, bv = b.d; else, bv = b; end
      if ga && ~gb && ~isscalar(bv) && ~isempty(bv)
        error('binary operator not implemented for ''ocl matrix'' by ''matrix'' operations (strict gpuArray stub)');
      elseif gb && ~ga && ~isscalar(av) && ~isempty(av)
        error('binary operator not implemented for ''matrix'' by ''ocl matrix'' operations (strict gpuArray stub)');
      end
      if islogical(av), av = double(av); end
      if islogical(bv), bv = double(bv); end
      x = op(av, bv);
    end
    function c = cellu(c)
      for k = 1:numel(c), if isa(c{k}, 'gpuArray'), c{k} = c{k}.d; end, end
    end
  end
end
