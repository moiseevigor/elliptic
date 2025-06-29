%Test function for theta_prime.m

% Test error handling - invalid j parameter
%!test
%! clear
%! try
%!     theta_prime(0, 0.5, 0.8); % j must be 1,2,3,4
%!     assert(false, "Invalid j parameter didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'J must be a scalar integer: 1, 2, 3, or 4.')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
%! clear
%! try
%!     theta_prime(5, 0.5, 0.8); % j must be 1,2,3,4
%!     assert(false, "Invalid j parameter didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'J must be a scalar integer: 1, 2, 3, or 4.')), ...
%!         'Unexpected error message: %s', err.message);
%! end

% Test error handling - module out of range
%!test
%! clear
%! try
%!     theta_prime(1, 0.5, 2); % module out of range
%!     assert(false, "Module out of range didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'M must be in the range 0 <= M <= 1')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
%! clear
%! try
%!     theta_prime(1, 0.5, -0.1); % module out of range
%!     assert(false, "Module out of range didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'M must be in the range 0 <= M <= 1')), ...
%!         'Unexpected error message: %s', err.message);
%! end

% Test error handling - complex input
%!test
%! clear
%! try
%!     theta_prime(1, 0.5i, 0.8); % complex input
%!     assert(false, "Complex input didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'Input arguments must be real')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
%! clear
%! try
%!     theta_prime(1, 0.5, 0.8i); % complex input
%!     assert(false, "Complex input didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'Input arguments must be real')), ...
%!         'Unexpected error message: %s', err.message);
%! end

% Test error handling - insufficient arguments
%!test
%! clear
%! try
%!     theta_prime(1, 0.5); % not enough arguments
%!     assert(false, "Insufficient arguments didn't throw an error.");
%! catch err
%!     % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'Not enough input arguments')), ...
%!         'Unexpected error message: %s', err.message);
%! end

% Test simple inputs - theta function values at z=0
%!test
%! clear
%! [th1, thp1] = theta_prime(1, 0, 0.5);
%! [th2, thp2] = theta_prime(2, 0, 0.5);
%! [th3, thp3] = theta_prime(3, 0, 0.5);
%! [th4, thp4] = theta_prime(4, 0, 0.5);
%! 
%! % At z=0: θ₁(0) = 0, θ₂(0) = 2q^(1/4), θ₃(0) = θ₃(0), θ₄(0) = θ₄(0)
%! assert(abs(th1) < 1e-12, 'θ₁(0) should be zero.');
%! assert(th2 > 0, 'θ₂(0) should be positive.');
%! assert(th3 > 0, 'θ₃(0) should be positive.');  
%! assert(th4 > 0, 'θ₄(0) should be positive.');

% Test specific values with known results
%!test
%! clear
%! j = 1; z = 0.5; m = 0.5;
%! [th, thp] = theta_prime(j, z, m);
%! 
%! % Verify using numerical differentiation
%! dz = 1e-8;
%! th_plus = theta(j, z + dz, m);
%! th_minus = theta(j, z - dz, m);
%! thp_numerical = (th_plus - th_minus) / (2 * dz);
%! 
%! % The derivative should match numerical differentiation within tolerance
%! assert(abs(thp - thp_numerical) < 1e-6, ...
%!     'Analytical derivative does not match numerical differentiation for θ₁.');

%!test
%! clear
%! j = 2; z = 0.3; m = 0.7;
%! [th, thp] = theta_prime(j, z, m);
%! 
%! % Verify using numerical differentiation
%! dz = 1e-8;
%! th_plus = theta(j, z + dz, m);
%! th_minus = theta(j, z - dz, m);
%! thp_numerical = (th_plus - th_minus) / (2 * dz);
%! 
%! % The derivative should match numerical differentiation within tolerance
%! assert(abs(thp - thp_numerical) < 1e-6, ...
%!     'Analytical derivative does not match numerical differentiation for θ₂.');

%!test
%! clear
%! j = 3; z = 0.2; m = 0.9;
%! [th, thp] = theta_prime(j, z, m);
%! 
%! % Verify using numerical differentiation
%! dz = 1e-8;
%! th_plus = theta(j, z + dz, m);
%! th_minus = theta(j, z - dz, m);
%! thp_numerical = (th_plus - th_minus) / (2 * dz);
%! 
%! % The derivative should match numerical differentiation within tolerance
%! assert(abs(thp - thp_numerical) < 1e-6, ...
%!     'Analytical derivative does not match numerical differentiation for θ₃.');

%!test
%! clear
%! j = 4; z = 0.4; m = 0.3;
%! [th, thp] = theta_prime(j, z, m);
%! 
%! % Verify using numerical differentiation
%! dz = 1e-8;
%! th_plus = theta(j, z + dz, m);
%! th_minus = theta(j, z - dz, m);
%! thp_numerical = (th_plus - th_minus) / (2 * dz);
%! 
%! % The derivative should match numerical differentiation within tolerance
%! assert(abs(thp - thp_numerical) < 1e-6, ...
%!     'Analytical derivative does not match numerical differentiation for θ₄.');

% Test array inputs with scalar j and m
%!test
%! clear
%! j = 1; 
%! z = [0.1, 0.2, 0.3, 0.4, 0.5];
%! m = 0.6;
%! [th, thp] = theta_prime(j, z, m);
%! 
%! assert(length(th) == length(z), 'Output size mismatch for array z input.');
%! assert(length(thp) == length(z), 'Output size mismatch for array z input.');
%! 
%! % Verify each element using numerical differentiation
%! for i = 1:length(z)
%!     dz = 1e-8;
%!     th_plus = theta(j, z(i) + dz, m);
%!     th_minus = theta(j, z(i) - dz, m);
%!     thp_numerical = (th_plus - th_minus) / (2 * dz);
%!     assert(abs(thp(i) - thp_numerical) < 1e-6, ...
%!         'Analytical derivative does not match numerical differentiation at z(%d).', i);
%! end

% Test array inputs with scalar j and z
%!test
%! clear
%! j = 2; 
%! z = 0.4;
%! m = [0.1, 0.3, 0.5, 0.7, 0.9];
%! [th, thp] = theta_prime(j, z, m);
%! 
%! assert(length(th) == length(m), 'Output size mismatch for array m input.');
%! assert(length(thp) == length(m), 'Output size mismatch for array m input.');
%! 
%! % Verify each element using numerical differentiation
%! for i = 1:length(m)
%!     dz = 1e-8;
%!     th_plus = theta(j, z + dz, m(i));
%!     th_minus = theta(j, z - dz, m(i));
%!     thp_numerical = (th_plus - th_minus) / (2 * dz);
%!     assert(abs(thp(i) - thp_numerical) < 1e-6, ...
%!         'Analytical derivative does not match numerical differentiation at m(%d).', i);
%! end

% Test array inputs with matching array sizes
%!test
%! clear
%! j = 3;
%! z = [0.1, 0.2, 0.3, 0.4];
%! m = [0.2, 0.4, 0.6, 0.8];
%! [th, thp] = theta_prime(j, z, m);
%! 
%! assert(length(th) == length(z), 'Output size mismatch for matched array inputs.');
%! assert(length(thp) == length(z), 'Output size mismatch for matched array inputs.');
%! 
%! % Verify each element using numerical differentiation
%! for i = 1:length(z)
%!     dz = 1e-8;
%!     th_plus = theta(j, z(i) + dz, m(i));
%!     th_minus = theta(j, z(i) - dz, m(i));
%!     thp_numerical = (th_plus - th_minus) / (2 * dz);
%!     assert(abs(thp(i) - thp_numerical) < 1e-6, ...
%!         'Analytical derivative does not match numerical differentiation at index %d.', i);
%! end

% Test limit cases - m approaching 0 and 1
%!test
%! clear
%! j = 1; z = 0.5;
%! 
%! % Test m approaching 0
%! m = 1e-10;
%! [th, thp] = theta_prime(j, z, m);
%! assert(isfinite(th) && isfinite(thp), 'Function should remain finite as m approaches 0.');
%! 
%! % Test m approaching 1
%! m = 1 - 1e-10;
%! [th, thp] = theta_prime(j, z, m);
%! assert(isfinite(th) && isfinite(thp), 'Function should remain finite as m approaches 1.');

% Test special values and symmetries 
%!test
%! clear
%! % Test θ₁ odd symmetry: θ₁(-z) = -θ₁(z)
%! j = 1; z = 0.3; m = 0.5;
%! [th_pos, thp_pos] = theta_prime(j, z, m);
%! [th_neg, thp_neg] = theta_prime(j, -z, m);
%! 
%! assert(abs(th_pos + th_neg) < 1e-12, 'θ₁ should satisfy odd symmetry: θ₁(-z) = -θ₁(z).');
%! assert(abs(thp_pos - thp_neg) < 1e-12, 'θ₁'' should satisfy even symmetry: θ₁''(-z) = θ₁''(z).');

% Test with tolerance parameter
%!test
%! clear
%! j = 1; z = 0.3; m = 0.7; tol = 1e-10;
%! [th, thp] = theta_prime(j, z, m, tol);
%! 
%! % Should produce same result as default tolerance for this simple case
%! [th_default, thp_default] = theta_prime(j, z, m);
%! assert(abs(th - th_default) < 1e-10, 'Results with custom tolerance should be close to default.');
%! assert(abs(thp - thp_default) < 1e-10, 'Results with custom tolerance should be close to default.');

% Benchmark test
%!test
%! clear
%! elapsedTime = [];
%! mem = [];
%! for i = 1:10
%!     z_vals = linspace(0, 1, 100);
%!     m_vals = linspace(0.1, 0.9, 100);
%!     [Z, M] = meshgrid(z_vals, m_vals);
%!     
%!     tic
%!     mem1 = whos();
%!     [th1, thp1] = theta_prime(1, Z, M);
%!     [th2, thp2] = theta_prime(2, Z, M);
%!     [th3, thp3] = theta_prime(3, Z, M);
%!     [th4, thp4] = theta_prime(4, Z, M);
%!     mem2 = whos();
%!     elapsedTime(i) = toc;
%!     mem(i) = sum([mem2.bytes]) - sum([mem1.bytes]);
%!     clear th1 th2 th3 th4 thp1 thp2 thp3 thp4 Z M;
%! end
% fprintf('\nAverage execution time for theta_prime calculations: %f seconds\n', mean(elapsedTime));
% fprintf('Average Mem: %f bytes\n', mean(mem));
%! assert(mean(elapsedTime) < 2.0, 'Average execution time for theta_prime calculations: %f seconds is greater than 2.0\n', mean(elapsedTime))
%! assert(mean(mem) < 1e7, 'Average memory used for theta_prime run: %f bytes is greater than 10MB\n', mean(mem)) 