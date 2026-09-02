function result = par_worker(func_name, varargin)
%PAR_WORKER  Generic parallel worker for Octave parcellfun.
%   Calls the named function with the given arguments and packs
%   multiple outputs into a cell array.
%
%   The parallel dispatch is switched OFF for the duration of the call and
%   restored afterwards.  Relying on parcellfun workers being fresh
%   processes (where the config defaults to parallel=false) was not
%   enough: whenever N is a multiple of chunk_size every chunk has exactly
%   chunk_size elements, so an in-process evaluation (a serial parcellfun
%   fallback, a shared-state worker, or plain testing) re-entered the
%   dispatch from inside the worker without bound and crashed Octave.

    was_parallel = elliptic_config('parallel');
    elliptic_config('parallel', false);
    unwind_protect
        result = par_worker_dispatch(func_name, varargin{:});
    unwind_protect_cleanup
        elliptic_config('parallel', was_parallel);
    end_unwind_protect
end


function result = par_worker_dispatch(func_name, varargin)
    switch func_name
        case 'elliptic12'
            [F, E, Z] = elliptic12(varargin{:});
            result = {F, E, Z};
        case 'elliptic3'
            result = elliptic3(varargin{:});
        case 'ellipj'
            [sn, cn, dn, am] = ellipj(varargin{:});
            result = {sn, cn, dn, am};
        case 'jacobiThetaEta'
            [Th, H] = jacobiThetaEta(varargin{:});
            result = {Th, H};
        case 'weierstrassP'
            result = weierstrassP(varargin{:});
        case 'weierstrassPPrime'
            result = weierstrassPPrime(varargin{:});
        case 'weierstrassZeta'
            result = weierstrassZeta(varargin{:});
        case 'weierstrassSigma'
            result = weierstrassSigma(varargin{:});
        case 'ellipticBD'
            [B, D, S] = ellipticBD(varargin{:});
            result = {B, D, S};
        case 'ellipticBDJ'
            if numel(varargin) >= 3
                [B, D, J] = ellipticBDJ(varargin{:});
                result = {B, D, J};
            else
                [B, D] = ellipticBDJ(varargin{:});
                result = {B, D};
            end
        case 'carlsonRF'
            result = carlsonRF(varargin{:});
        case 'carlsonRD'
            result = carlsonRD(varargin{:});
        case 'carlsonRJ'
            result = carlsonRJ(varargin{:});
        case 'carlsonRC'
            result = carlsonRC(varargin{:});
        case 'jacobiEDJ'
            if numel(varargin) >= 3
                [Eu, Du, Ju] = jacobiEDJ(varargin{:});
                result = {Eu, Du, Ju};
            else
                [Eu, Du] = jacobiEDJ(varargin{:});
                result = {Eu, Du};
            end
        case 'cel'
            result = cel(varargin{:});
        otherwise
            error('par_worker: unknown function %s', func_name);
    end
end
