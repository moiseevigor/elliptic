function [C, ia, ic] = uniquetol_compat(A, tol)
% UNIQUETOL_COMPAT  Compatibility wrapper for uniquetol.
%   [C, IA, IC] = UNIQUETOL_COMPAT(A, TOL) returns unique values within
%   tolerance TOL. Uses MATLAB's built-in uniquetol when available,
%   otherwise falls back to a sort-based implementation for Octave.

if exist('uniquetol', 'builtin') || exist('uniquetol', 'file')
    [C, ia, ic] = uniquetol(A, tol);
else
    % Octave-compatible fallback using sort for O(n log n) performance
    A = A(:).';
    n = length(A);
    [sortedA, sortIdx] = sort(A);

    % Group sorted elements: compare each to its group representative
    group = 1;
    groupOf = ones(1, n);       % group assignment in sorted order
    repVal = sortedA(1);        % representative value for current group
    ia_sorted = [1];            % indices in sorted array of group representatives

    for i = 2:n
        if abs(sortedA(i) - repVal) > tol * max(1, abs(repVal))
            group = group + 1;
            repVal = sortedA(i);
            ia_sorted(end+1) = i;
        end
        groupOf(i) = group;
    end

    % Map back to original order
    ic = zeros(1, n);
    ic(sortIdx) = groupOf;

    ia = sort(sortIdx(ia_sorted));
    C = A(ia);
end
