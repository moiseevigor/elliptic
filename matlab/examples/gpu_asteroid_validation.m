function gpu_asteroid_validation()
% GPU_ASTEROID_VALIDATION  Summary of Keplerian-vs-Horizons error samples.
%
% Reads examples/gpu-asteroid-swarm/data/validation_sample.json produced by
% scripts/validate_against_horizons.py and prints a table.

  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  jsonPath = fullfile(root, 'examples', 'gpu-asteroid-swarm', 'data', 'validation_sample.json');
  if ~exist(jsonPath, 'file')
    error('Run scripts/validate_against_horizons.py first (missing %s)', jsonPath);
  end
  payload = jsondecode(fileread(jsonPath));
  fprintf('Validation date: %s  (JD %.1f)\n', payload.date, payload.jd);
  fprintf('%-26s %9s %12s\n', 'Object', 'dT [yr]', '|dr| [AU]');
  fprintf('%s\n', repmat('-', 1, 50));
  for k = 1:numel(payload.samples)
    s = payload.samples(k);
    fprintf('%-26s %+9.2f %12.4e\n', s.name, s.dt_years, s.err_au);
  end
end
