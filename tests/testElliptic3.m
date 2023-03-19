%Test function for elliptic3.m

%!test
% Test error handling
%! try
%!     elliptic3(0, 2, 1); % module out of range
%!     assert(false, "Module out of range didn't throw an error.");
%! catch err
    % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'M and C must be in the range [0, 1].')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
%! try
%!     elliptic3(0, 0.5i, 0.5); % complex input
%!     assert(false, "Complex input didn't throw an error.");
%! catch err
    % Verify that the error message contains the expected string
%!     assert(~isempty(strfind(err.message, 'Input arguments must be real.')), ...
%!         'Unexpected error message: %s', err.message);
%! end

%!test
% Test some simple inputs
%! Pi = elliptic3(0, 0.5, 1);
%! assert(abs(Pi - 0) < 1e-12, 'Pi value is incorrect.');

%!test
% Test the output of elliptic3 for some inputs
%! Pi = elliptic3(pi/4, 0.5, 0.5);
%! assert(abs(Pi - 0.919022739165697) < 1e-10, 'Unexpected value for K');
    
%!test
% Test a range of inputs
%! [phi,alpha,c] = meshgrid(0:25:90, 0:25:90, 0:0.4:1); 
%! Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
%! expectedPi = [
%! [
%!         0   0.436332312998582   0.872664625997165   1.308996938995747
%!         0   0.438747923080541   0.890543879394498   1.360834669550633
%!         0   0.444551505567130   0.940075683068702   1.534546187668256
%!         0   0.449816440114804   0.997105354291443   1.871453962422796
%!     ]
%!     [
%!         0   0.447481674284817   0.962368261625533   1.597941344935833
%!         0   0.450008007383737   0.983487876637982   1.668946178349462
%!         0   0.456080142788474   1.042322814335382   1.911003257311173
%!         0   0.461591818685859   1.110665166734768   2.397752070491359
%!     ]
%!     [
%!         0   0.459719055660277   1.094942500776336   2.305386400623875
%!         0   0.462369007025916   1.121161300350779   2.428736708263183
%!         0   0.468741143521980   1.194760969067305   2.863258220259662
%!         0   0.474528668429618   1.281307083403909   3.803705851444642
%!     ]
%! ];