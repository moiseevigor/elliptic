function r = sub_kpi(u, k)
%SUB_KPI  r = u - k*pi with the product formed exactly (three-term split of pi).
%   R = SUB_KPI(U, K) returns U - K*PI for integer-valued K, accurate to
%   eps*|R| rather than eps*|U|.  PI is split as PI_A + PI_B + PI_C where
%   PI_A and PI_B carry 25 significant bits each, so K*PI_A and K*PI_B are
%   exact in double for |K| < 2^28 (|U| < 8e8); the remaining K*PI_C rounds
%   at eps*|K|*1.6e-8, far below eps*|R|.  Using double(pi) as the leading
%   term (the previous "Cody-Waite" split) does not help: K*double(pi)
%   already rounds by eps*|U| (2.3e-10 at U = 1e6), which Jacobi Zeta and E
%   inherit.
%
%   PI_A = 0x1.921fb5p+1, PI_B = 0x1.110b46p-26, PI_C = pi - PI_A - PI_B
%   (residual 1.3e-24).  Works elementwise on host or GPU arrays.
r = ((u - k .* 3.1415926218032837) - k .* 1.5893254712295857e-08) - k .* 1.5893254834760535e-08;
end
