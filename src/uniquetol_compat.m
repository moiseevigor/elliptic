function [C, ia, ic] = uniquetol_compat(A, tol)
% UNIQUETOL_COMPAT  Compatibility wrapper for uniquetol.
%   [C, IA, IC] = UNIQUETOL_COMPAT(A, TOL) returns unique values within
%   tolerance TOL. Uses MATLAB's built-in uniquetol when available,
%   otherwise falls back to a vectorized sort-based implementation for Octave.

if exist('uniquetol', 'builtin') || exist('uniquetol', 'file')
    [C, ia, ic] = uniquetol(A, tol);
else
    % Vectorized Octave-compatible fallback
    A = A(:).';
    [sortedA, sortIdx] = sort(A);

    % Mark the start of each new group: first element is always a new group,
    % then any element that differs from its predecessor by more than tolerance.
    diffs = [true, abs(diff(sortedA)) > tol * max(1, abs(sortedA(1:end-1)))];
    groupOf = cumsum(diffs);   % group index for each sorted element (1-based)

    % Map group indices back to original order
    ic = zeros(1, length(A));
    ic(sortIdx) = groupOf;

    % Representative for each group: first occurrence in sorted order
    ia_sorted = find(diffs);   % positions in sorted array where each group starts
    ia = sortIdx(ia_sorted);   % corresponding positions in original array
    C = sortedA(ia_sorted);    % unique values in sorted ascending order
end
