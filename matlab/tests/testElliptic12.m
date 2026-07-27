%Test function for elliptic12.m

% Test error handling
%!test
%! clear
%! try
%!     elliptic12(0, 2); % module out of range
%!     assert(false, "Module out of range didn't throw an error.");
%! catch err
    % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'M must be in the range 0 <= M <= 1')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
%! clear
%! try
%!     elliptic12(0, 0.5i); % complex input
%!     assert(false, "Complex input didn't throw an error.");
%! catch err
    % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'Input arguments must be real. Use ELLIPTIC12i for complex arguments')), ...
%!         'Unexpected error message: %s', err.message);
%! end

% Test some simple inputs
%!test
%! clear
%! [F,E,Z] = elliptic12(0, 0.5);
%! assert(abs(F - 0) < 1e-12, 'F value is incorrect.');
%! assert(abs(E - 0) < 1e-12, 'E value is incorrect.');
%! assert(abs(Z - 0) < 1e-12, 'Z value is incorrect.');

% Test the output of elliptic12 for some inputs
%!test
%! clear
%! [F,E,Z] = elliptic12(1000*pi/e, 0.5);
%! assert(abs(F - 1364.215673994739) < 1e-10, 'Unexpected value for F');
%! assert(abs(E - 993.6995958059659) < 1e-10, 'Unexpected value for E');
%! assert(abs(Z - (-9.508521098575250e-02)) < 1e-10, 'Unexpected value for Z');
    
% Test a range of inputs
%!test
%! clear
%! [phi,alpha] = meshgrid(0:15:90, 0:20:90); 
%! [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
%! expectedF = [
%!     0   0.261799387799149   0.523598775598299   0.785398163397448   1.047197551196598   1.308996938995747   1.570796326794897
%!     0   0.262145681692449   0.526283990562203   0.793981429961732   1.065968913708522   1.341839009622797   1.620025899124204
%!     0   0.263033690353369   0.533427451037688   0.818147652199543   1.122556669683242   1.447669376453414   1.786769134885021
%!     0   0.264063548276829   0.542229109803553   0.851223749071185   1.212596615254979   1.649178665655556   2.156515647499643
%!     0   0.264747663195422   0.548425344542772   0.877408330405700   1.301353213761152   1.946822305295344   3.153385251887838
%! ];
%! assert(norm(F-expectedF) < 1e-12, 'F value is incorrect.')
%! expectedE = [
%!     0   0.261799387799149   0.523598775598299   0.785398163397448   1.047197551196598   1.308996938995747   1.570796326794897
%!     0   0.261453912906691   0.520937696463332   0.776974018512720   1.028972213953050   1.277421529463950   1.523799205259775
%!     0   0.260575450957998   0.514088617513829   0.754888085407355   0.980134299660981   1.191010358808463   1.393140248523812
%!     0   0.259569955098381   0.506092072465726   0.728224155457347   0.918393294316326   1.075856687976766   1.211056027568461
%!     0   0.258909826877359   0.500742319367686   0.709723805114461   0.872755203912916   0.981407813910667   1.040114395706011
%! ];
%! assert(norm(E-expectedE) < 1e-12, 'E value is incorrect.')
%! expectedZ = [
%!     0                   0                   0                   0                   0                   0                   0
%!     0   0.014879224414869   0.025914050857933   0.030153691368092   0.026319982023431   0.015285400927267   0.000000000000000
%!     0   0.055488619316143   0.098176322411480   0.116980030437629   0.104879154913245   0.062265499992439   0.000000000000000
%!     0   0.111277151300479   0.201587056463619   0.250193934198656   0.237423303852160   0.149710955672358   0.000000000000000
%!     0   0.171585306172434   0.319849389939710   0.420318939400386   0.443516115310938   0.339266830849054   0.000000000000000
%! ];
%! assert(norm(Z-expectedZ) < 1e-12, 'Z value is incorrect.')

% Issue #35: quasi-period past pi/2 must be 2*k*K(m) / 2*k*E(m), not k*K / k*E
%!test
%! clear
%! m = 0.4;
%! [K,Ec] = ellipke(m);
%! phi = [-1.2 1.2 pi/2 1.6 2.7 pi 4.0 3*pi/2 7.0];
%! [F,E] = elliptic12(phi, m*ones(size(phi)));
%! for k = 1:numel(phi)
%!     Fq = quadgk(@(t) 1./sqrt(1-m*sin(t).^2), 0, phi(k), 'AbsTol', 1e-13);
%!     Eq = quadgk(@(t)  sqrt(1-m*sin(t).^2),   0, phi(k), 'AbsTol', 1e-13);
%!     assert(abs(F(k)-Fq) < 1e-11, 'F(%g) = %.12f, expected %.12f', phi(k), F(k), Fq);
%!     assert(abs(E(k)-Eq) < 1e-11, 'E(%g) = %.12f, expected %.12f', phi(k), E(k), Eq);
%! end
%! assert(abs(elliptic12(pi/2 + pi, m) - (elliptic12(pi/2, m) + 2*K)) < 1e-12, 'F is not 2K-quasi-periodic.')

% Issue #35: elliptic12i must stay continuous across Re(u) = pi/2 below the
% branch point at pi/2 + i*acosh(1/sqrt(m)), and match elliptic12 on the real axis
%!test
%! clear
%! m = 0.4;
%! psi = 0.5;                                  % < acosh(1/sqrt(0.4)) = 1.0317
%! phi = linspace(0.2, 3.0, 141);
%! [Fi,Ei] = elliptic12i(phi + 1i*psi, m*ones(size(phi)));
%! step = max(abs(diff(Fi)));                  % ~0.03 when smooth, ~K when branching
%! assert(step < 0.05, 'F jumps by %g across Re(u) = pi/2.', step);
%! step = max(abs(diff(Ei)));
%! assert(step < 0.05, 'E jumps by %g across Re(u) = pi/2.', step);
%! phi = [0.3 1.2 pi/2 2.0 3.4];               % real u: must agree with elliptic12
%! [Fr,Er] = elliptic12(phi, m*ones(size(phi)));
%! [Fc,Ec] = elliptic12i(phi, m*ones(size(phi)));
%! assert(norm(Fc-Fr) < 1e-9 && norm(Ec-Er) < 1e-9, 'elliptic12i disagrees with elliptic12 on the real axis.')

% Issue #35: at Re(u) = pi/2 exactly the imaginary part used to collapse to zero.
% Closed form, valid for |psi| < acosh(1/sqrt(m)):
%   F(pi/2 + i*psi | m) = K(m) + i*F(mu | 1-m),  sin(mu) = tanh(psi)/sqrt(1-m)
%!test
%! clear
%! for m = [0.1 0.4 0.75]
%!     for psi = [0.01 0.05 0.3]
%!         mu   = asin(tanh(psi)/sqrt(1-m));
%!         want = ellipke(m) + 1i*elliptic12(mu, 1-m);
%!         got  = elliptic12i(pi/2 + 1i*psi, m);
%!         assert(abs(got-want) < 1e-12, ...
%!             'F(pi/2+%gi | %g) = %.12f%+.12fi, expected %.12f%+.12fi', ...
%!             psi, m, real(got), imag(got), real(want), imag(want));
%!     end
%! end

% Benchmark time and memory
%!test
%! clear
%! elapsedTime = [];
%! mem = [];
%! for i=1:10
%!     [phi,alpha] = meshgrid(0:0.5:720, 0:0.5:90); 
%!     tic
%!     mem1 = whos();
%!     [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
%!     mem2 = whos();
%!     elapsedTime(i) = toc;
%!     mem(i) = sum([mem2.bytes]) - sum([mem1.bytes]);
%!     clear F E Z phi alpha;
%! end
% fprintf('\nAverage execution time for elliptic12 calculations: %f seconds\n', mean(elapsedTime));
% fprintf('Average Mem: %f\n', mean(mem));
%! thresh = 0.6;  % Octave uses vectorized uniquetol_compat (~0.16s); MATLAB built-in uniquetol is faster (~0.05s). CI machines can be 2x slower than dev.
%! assert(mean(elapsedTime) < thresh, 'Average execution time for elliptic12 calculations: %f seconds is greater than %f\n', mean(elapsedTime), thresh)
%! assert(mean(mem) < 6259742.2, 'Average memory used for elliptic12 run: %f bytes is greater than 6259742.1\n', mean(mem))

