function result = par_worker(func_name, varargin)
%PAR_WORKER  Generic parallel worker for Octave parcellfun.
%   Calls the named function with the given arguments and packs
%   multiple outputs into a cell array. Workers run in fresh
%   processes where elliptic_config defaults to parallel=false,
%   preventing recursive parallelism.

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
        otherwise
            error('par_worker: unknown function %s', func_name);
    end
