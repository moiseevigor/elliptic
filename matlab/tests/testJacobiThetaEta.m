% Tests for jacobiThetaEta.m - Jacobi Theta and Eta functions

%!test
%! % At u=0: H(0,m) = 0 for any m, Th(0,m) should be nonzero
%! [Th, H] = jacobiThetaEta(0, 0.5);
%! assert(abs(H) < 1e-14, 'H(0, 0.5) should be zero.');
%! assert(Th > 0, 'Th(0, 0.5) should be positive.');

%!test
%! % Special case m=0: Th should be 1, H should be sqrt(sqrt(m))*sin(u) ~ 0
%! [Th, H] = jacobiThetaEta(0.5, 0);
%! assert(abs(Th - 1) < 1e-14, 'Th(u, 0) should be 1.');
%! assert(abs(H) < 1e-14, 'H(u, 0) should be ~0 for m=0.');

%!test
%! % Special case m=1: should return NaN
%! [Th, H] = jacobiThetaEta(0.5, 1);
%! assert(isnan(Th), 'Th(u, 1) should be NaN.');
%! assert(isnan(H), 'H(u, 1) should be NaN.');

%!test
%! % Vector input
%! u = [0 0.3 0.6 0.9];
%! m = [0.5 0.5 0.5 0.5];
%! [Th, H] = jacobiThetaEta(u, m);
%! assert(length(Th) == 4, 'Output should have same size as input.');
%! assert(abs(H(1)) < 1e-14, 'H(0, m) should be zero.');

%!test
%! % Mixed special cases: m=0 and m=1 in same call
%! u = [0.5 0.5 0.5];
%! m = [0   0.5 1  ];
%! [Th, H] = jacobiThetaEta(u, m);
%! assert(abs(Th(1) - 1) < 1e-14, 'Th(u, 0) should be 1.');
%! assert(isnan(Th(3)), 'Th(u, 1) should be NaN.');
%! assert(isnan(H(3)), 'H(u, 1) should be NaN.');

%!test
%! % Error on complex input
%! try
%!     jacobiThetaEta(1+1i, 0.5);
%!     assert(false, 'Should have thrown error for complex input.');
%! catch err
%!     assert(~isempty(strfind(err.message, 'real')), ...
%!         'Unexpected error: %s', err.message);
%! end

%!test
%! % Error on m out of range
%! try
%!     jacobiThetaEta(0.5, 1.5);
%!     assert(false, 'Should have thrown error for m > 1.');
%! catch err
%!     assert(~isempty(strfind(err.message, 'range')), ...
%!         'Unexpected error: %s', err.message);
%! end
