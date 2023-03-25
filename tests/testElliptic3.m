%Test function for elliptic3.m

% Test error handling
%!test
%! clear
%! try
%!     elliptic3(0, 2, 1); % module out of range
%!     assert(false, "Module out of range didn't throw an error.");
%! catch err
    % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'M and C must be in the range [0, 1].')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
%! clear
%! try
%!     elliptic3(0, 0.5i, 0.5); % complex input
%!     assert(false, "Complex input didn't throw an error.");
%! catch err
    % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'Input arguments must be real.')), ...
%!         'Unexpected error message: %s', err.message);
%! end

% Test some simple inputs
%!test
%! clear
%! Pi = elliptic3(0, 0.5, 1);
%! assert(abs(Pi - 0) < 1e-12, 'Pi value is incorrect.');

% Test the output of elliptic3 for some inputs
%!test
%! clear
%! Pi = elliptic3(pi/4, 0.5, 0.5);
%! assert(abs(Pi - 0.919022739165697) < 1e-10, 'Unexpected value for K');
    
% Test a range of inputs
%!test
%! clear
%! [phi,alpha,c] = meshgrid(0:25:90, 0:25:90, 0:0.4:1); 
%! Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
%! expectedPi = [
%!     0                   0                   0
%!     0                   0                   0
%!     0                   0                   0
%!     0                   0                   0
%!     0.436332312998582   0.447481674284817   0.459719055660277
%!     0.438747923080541   0.450008007383737   0.462369007025916
%!     0.444551505567130   0.456080142788474   0.468741143521980
%!     0.449816440114804   0.461591818685859   0.474528668429618
%!     0.872664625997165   0.962368261625533   1.094942500776336
%!     0.890543879394498   0.983487876637982   1.121161300350779
%!     0.940075683068702   1.042322814335382   1.194760969067305
%!     0.997105354291443   1.110665166734768   1.281307083403909
%!     1.308996938995747   1.597941344935833   2.305386400623875
%!     1.360834669550633   1.668946178349462   2.428736708263183
%!     1.534546187668256   1.911003257311173   2.863258220259662
%!     1.871453962422796   2.397752070491359   3.803705851444642
%! ];
%! assert(size(Pi) == [4 4 3], 'Pi size is incorrect.')
%! Pi = reshape(Pi, [], size(Pi, 3));
%! assert(norm(Pi-expectedPi) < 1e-12, 'Pi value is incorrect.')

% Benchmark time and memory
%!test
%! clear
%! elapsedTime = [];
%! mem = [];
%! for i=1:10
%!     [phi,alpha,c] = meshgrid(0:0.5:90, 0:0.5:90, 0:0.1:1); 
%!     tic
%!     mem1 = whos();
%!     Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
%!     mem2 = whos();
%!     elapsedTime(i) = toc;
%!     mem(i) = sum([mem2.bytes]) - sum([mem1.bytes]);
%!     clear Pi phi alpha;
%! end
% fprintf('\nAverage execution time for elliptic3 calculations: %f seconds\n', mean(elapsedTime));
% fprintf('Average Mem: %f\n', mean(mem));
%! assert(mean(elapsedTime) < 0.15, 'Average execution time for elliptic3 calculations: %f seconds is greater than 0.15\n', mean(elapsedTime))
%! assert(mean(mem) < 2883013.7, 'Average memory used for elliptic3 run: %f bytes is greater than 2883013.7\n', mean(mem))

