% Tests for agm.m - Arithmetic Geometric Mean

%!test
%! % AGM(1, sqrt(2)/2) should converge to known value
%! % AGM(1, 1/sqrt(2)) = pi / (2 * K(1/2)) where K is complete elliptic integral
%! [a,b,c,n] = agm(1, sqrt(2)/2, 1 - sqrt(2)/2);
%! mn = double(max(n));
%! assert(abs(a(mn+1) - b(mn+1)) < 1e-15, 'AGM did not converge: a and b should be equal.');

%!test
%! % AGM(a,a) = a for any a
%! [a,b,c,n] = agm(5, 5, 0);
%! assert(abs(a(1) - 5) < 1e-15, 'AGM(a,a) should equal a.');
%! assert(abs(c(1)) < 1e-15, 'C should be zero when a==b.');

%!test
%! % Vector input: AGM of multiple pairs
%! a0 = [1 1 1];
%! b0 = [0.5 0.7 0.9];
%! c0 = a0 - b0;
%! [a,b,c,n] = agm(a0, b0, c0);
%! for k = 1:3
%!     mn = double(n(k));
%!     assert(abs(c(mn+1,k)) < 1e-14, sprintf('AGM did not converge for pair %d.', k));
%! end

%!test
%! % AGM(1,0) should converge to 0
%! [a,b,c,n] = agm(1, 0, 1);
%! mn = double(max(n));
%! assert(abs(b(mn+1)) < 1e-15, 'AGM(1,0) should converge to 0.');
