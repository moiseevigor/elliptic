% generate_data.m — Pre-compute reference data for the interactive page.
% Run from the repo root: octave --no-gui examples/arclength-celestial-mechanics/scripts/generate_data.m
% Writes JSON to examples/arclength-celestial-mechanics/data/

addpath('src');

OUT = 'examples/arclength-celestial-mechanics/data';

% ── 1. E(phi|m) reference table ──────────────────────────────────────────────
phi_v = linspace(0, pi/2, 100);
m_v   = [0.0, 0.2, 0.4, 0.6, 0.8, 0.95];

fid = fopen(fullfile(OUT, 'ellipse_arcs.json'), 'w');
fprintf(fid, '{\n');
fprintf(fid, '  "phi": ['); fprintf(fid, '%.10g,', phi_v(1:end-1)); fprintf(fid, '%.10g],\n', phi_v(end));
fprintf(fid, '  "m_values": ['); fprintf(fid, '%.4g,', m_v(1:end-1)); fprintf(fid, '%.4g],\n', m_v(end));
fprintf(fid, '  "E": [\n');
for mi = 1:length(m_v)
  [~, E_v] = elliptic12(phi_v, m_v(mi).*ones(size(phi_v)));
  fprintf(fid, '    ['); fprintf(fid, '%.15g,', E_v(1:end-1)); fprintf(fid, '%.15g]', E_v(end));
  if mi < length(m_v), fprintf(fid, ','); end
  fprintf(fid, '\n');
end
fprintf(fid, '  ],\n');
fprintf(fid, '  "F": [\n');
for mi = 1:length(m_v)
  [F_v] = elliptic12(phi_v, m_v(mi).*ones(size(phi_v)));
  fprintf(fid, '    ['); fprintf(fid, '%.15g,', F_v(1:end-1)); fprintf(fid, '%.15g]', F_v(end));
  if mi < length(m_v), fprintf(fid, ','); end
  fprintf(fid, '\n');
end
fprintf(fid, '  ]\n');
fprintf(fid, '}\n');
fclose(fid);

% ── 2. Kepler orbit data ──────────────────────────────────────────────────────
e_v   = [0.0, 0.3, 0.5, 0.7, 0.9];
theta = linspace(0, 2*pi, 361);

fid = fopen(fullfile(OUT, 'kepler_orbits.json'), 'w');
fprintf(fid, '{\n');
fprintf(fid, '  "e_values": ['); fprintf(fid, '%.2g,', e_v(1:end-1)); fprintf(fid, '%.2g],\n', e_v(end));
fprintf(fid, '  "theta": ['); fprintf(fid, '%.8g,', theta(1:end-1)); fprintf(fid, '%.8g],\n', theta(end));
fprintf(fid, '  "r_normalized": [\n');
a = 1;
for ei = 1:length(e_v)
  e = e_v(ei);
  r = a*(1-e^2)./(1 + e.*cos(theta));
  fprintf(fid, '    ['); fprintf(fid, '%.10g,', r(1:end-1)); fprintf(fid, '%.10g]', r(end));
  if ei < length(e_v), fprintf(fid, ','); end
  fprintf(fid, '\n');
end
fprintf(fid, '  ],\n');

% Arc length from perihelion via eccentric anomaly
fprintf(fid, '  "arc_L_normalized": [\n');
for ei = 1:length(e_v)
  e = e_v(ei);
  m = e^2;
  % Eccentric anomaly from true anomaly
  E_anom = 2.*atan(sqrt((1-e)./(1+e)).*tan(theta./2));
  E_anom(theta > pi) = E_anom(theta > pi) + 2*pi;
  % Arc length: a*[E_complete - E(pi/2 - E_anom, m)]
  % Use the quarter-period approach
  [~, Ec] = ellipke(m);
  Nq = floor(E_anom./(pi/2));
  r_frac = E_anom - Nq.*(pi/2);
  L = Nq .* Ec;
  for k = 0:max(Nq)
    idx_e = (Nq == k) & (mod(k,2)==0);
    idx_o = (Nq == k) & (mod(k,2)==1);
    if any(idx_e)
      [~, Ev] = elliptic12(pi/2 - r_frac(idx_e), m.*ones(1,sum(idx_e)));
      L(idx_e) = L(idx_e) + Ec - Ev;
    end
    if any(idx_o)
      [~, Ev] = elliptic12(r_frac(idx_o), m.*ones(1,sum(idx_o)));
      L(idx_o) = L(idx_o) + Ev;
    end
  end
  fprintf(fid, '    ['); fprintf(fid, '%.10g,', L(1:end-1)); fprintf(fid, '%.10g]', L(end));
  if ei < length(e_v), fprintf(fid, ','); end
  fprintf(fid, '\n');
end
fprintf(fid, '  ]\n');
fprintf(fid, '}\n');
fclose(fid);

% ── 3. Validation table: F, E, B, D at specific (phi, m) ─────────────────────
phi_tab = [pi/6, pi/4, pi/3, pi/2, pi/6, pi/4, pi/3, pi/2];
m_tab   = [0.3,  0.3,  0.3,  0.3,  0.7,  0.7,  0.7,  0.7 ];
[F_tab, E_tab] = elliptic12(phi_tab, m_tab);
[B_tab, D_tab] = ellipticBDJ(phi_tab, m_tab);

fid = fopen(fullfile(OUT, 'validation.json'), 'w');
fprintf(fid, '[\n');
for i = 1:length(phi_tab)
  fprintf(fid, '  {"phi": %.10g, "m": %.1f, "F": %.15g, "E": %.15g, "B": %.15g, "D": %.15g}', ...
    phi_tab(i), m_tab(i), F_tab(i), E_tab(i), B_tab(i), D_tab(i));
  if i < length(phi_tab), fprintf(fid, ','); end
  fprintf(fid, '\n');
end
fprintf(fid, ']\n');
fclose(fid);

fprintf('Data written to %s\n', OUT);
