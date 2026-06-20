function test_suite = test_qrupdate_gs %#ok<*STOUT>
    %
    %   Unit tests for the qrupdate_gs utility.
    %
    %   Modified from:
    %   https://github.com/Remi-Gau/template_matlab_analysis
    %
    %   Copyright:
    %   ----------
    %
    %   This file is part of the KrySBAS MATLAB Toolbox.
    %
    %   Copyright 2023 CC&MA - NIDTec - FP - UNA
    %
    %   KrySBAS is free software: you can redistribute it and/or modify it under
    %   the terms of the GNU General Public License as published by the Free
    %   Software Foundation, either version 3 of the License, or (at your
    %   option) any later version.
    %
    %   KrySBAS is distributed in the hope that it will be useful, but WITHOUT
    %   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    %   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    %   for more details.
    %
    %   You should have received a copy of the GNU General Public License along
    %   with this file.  If not, see <http://www.gnu.org/licenses/>.
    %

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end

    initTestSuite;

end

% =========================================================================
% ----> Full one-shot factorisation
% =========================================================================

function test_full_factorisation_correctness()
    % qrupdate_gs called all at once (empty initial Q, R) must reproduce
    % the standard thin QR factorisation: A = Q*R with Q'*Q = I and
    % R upper triangular.

    rng(0);
    A = randn(10, 4);

    [q_out, r_out] = qrupdate_gs(A, [], []);

    % A = Q * R
    assertElementsAlmostEqual(q_out * r_out, A);

    % Q has orthonormal columns
    assertElementsAlmostEqual(q_out' * q_out, eye(4));

    % R is upper triangular
    assertEqual(tril(r_out, -1), zeros(4));
end

function test_full_factorisation_diagonal_positive()
    % The diagonal entries of R from qrupdate_gs are always non-negative
    % (we divide by norm(w), which is positive).

    rng(1);
    A = randn(8, 3);

    [~, r_out] = qrupdate_gs(A, [], []);

    assert(all(diag(r_out) >= 0));
end

% =========================================================================
% ----> Incremental column-by-column update
% =========================================================================

function test_incremental_matches_full()
    % Building the QR column by column must give the same result as
    % calling qrupdate_gs on the full matrix at once.

    rng(2);
    A = randn(10, 5);

    % One-shot reference
    [q_ref, r_ref] = qrupdate_gs(A, [], []);

    % Column-by-column incremental build
    q_inc = [];
    r_inc = [];
    for j = 1:5
        [q_inc, r_inc] = qrupdate_gs(A(:, 1:j), q_inc, r_inc);
    end

    assertElementsAlmostEqual(q_inc, q_ref);
    assertElementsAlmostEqual(r_inc, r_ref);
end

function test_incremental_orthogonality_preserved()
    % After each incremental update the leading columns of Q must remain
    % orthonormal, i.e. Q(:,1:j)' * Q(:,1:j) = I_j at every step j.

    rng(3);
    A = randn(12, 6);

    q_inc = [];
    r_inc = [];
    for j = 1:6
        [q_inc, r_inc] = qrupdate_gs(A(:, 1:j), q_inc, r_inc);
        assertElementsAlmostEqual(q_inc' * q_inc, eye(j));
    end
end

% =========================================================================
% ----> Growing row count (the GMRES-DR use case)
% =========================================================================

function test_row_growth_between_calls()
    % In gmres_dr the Hessenberg gains one row at every Arnoldi step:
    % at step j the matrix is (j+1)-by-j.  This test mimics that pattern
    % and verifies that qrupdate_gs handles the row-count increase correctly.
    %
    % Concretely we check that A(1:j+1, 1:j) = Q * R at each step j,
    % using a fixed (tall) matrix A whose leading submatrices are passed
    % one column at a time.

    rng(4);
    m = 8;  % maximum number of columns
    A = randn(m + 1, m);  % (m+1)-by-m tall matrix

    q_inc = [];
    r_inc = [];
    for j = 1:m
        % Pass the (j+1)-by-j leading submatrix: one more row than last time
        a_sub = A(1:j + 1, 1:j);
        [q_inc, r_inc] = qrupdate_gs(a_sub, q_inc, r_inc);

        % Verify the factorisation is correct at each step
        assertElementsAlmostEqual(q_inc * r_inc, a_sub);
        assertElementsAlmostEqual(q_inc' * q_inc, eye(j));
    end
end

function test_row_growth_residual_estimate()
    % The key use of qrupdate_gs in gmres_dr is to compute the
    % least-squares residual norm cheaply.  This test checks that
    %
    %   d     = R \ (Q' * c)
    %   res   = norm(c - H * d)
    %
    % with an incrementally maintained QR matches the result from a
    % direct pinv solve on the full system.

    rng(5);
    m = 6;
    H = randn(m + 1, m);  % simulated Hessenberg
    c = randn(m + 1, 1);  % simulated RHS (Vr in GMRES-DR)

    % Incremental QR
    q_inc = [];
    r_inc = [];
    for j = 1:m
        h_sub = H(1:j + 1, 1:j);
        [q_inc, r_inc] = qrupdate_gs(h_sub, q_inc, r_inc);
    end

    d_inc = r_inc \ (q_inc' * c);
    res_inc = norm(c - H * d_inc);

    % Reference: direct pinv
    d_ref = pinv(H) * c;
    res_ref = norm(c - H * d_ref);

    assertElementsAlmostEqual(res_inc, res_ref);
end

% =========================================================================
% ----> Edge cases
% =========================================================================

function test_single_column()
    % A single-column tall matrix should produce a unit-vector Q and a
    % scalar R equal to the column norm.

    a = [3; 4];  % 2-by-1

    [q_out, r_out] = qrupdate_gs(a, [], []);

    assertElementsAlmostEqual(r_out, 5);
    assertElementsAlmostEqual(q_out, [3; 4] / 5);
    assertElementsAlmostEqual(q_out' * q_out, 1);
end
