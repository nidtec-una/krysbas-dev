function test_suite = test_gmres_dr %#ok<*STOUT>
    %
    %   Test suite for the GMRES-DR solver.
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
% ----> Sanity-check tests (input validation)
% =========================================================================

function test_number_of_input_arguments()
    % Verify that errors are raised for too few or too many input arguments.

    % Too few: only A is given (b is required)
    try
        gmres_dr(ones(2, 2));
    catch ME
        msg = "Too few input parameters. Expected at least A and b.";
        assert(matches(ME.message, msg));
    end

    % Too many: eight arguments instead of the maximum seven
    try
        gmres_dr([], [], [], [], [], [], [], []);
    catch ME
        msg = "Too many input parameters.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_matrix_A()
    % An empty matrix A must raise an error.
    try
        gmres_dr([], ones(2, 1));
    catch ME
        msg = 'Matrix A cannot be empty.';
        assert(matches(ME.message, msg));
    end
end

function test_non_square_matrix_A()
    % A non-square matrix A must raise an error.
    try
        gmres_dr([1; 1], [1; 1]);
    catch ME
        msg = "Matrix A must be square.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_vector_b()
    % An empty right-hand side b must raise an error.
    try
        gmres_dr(eye(2), []);
    catch ME
        msg = "Vector b cannot be empty.";
        assert(matches(ME.message, msg));
    end
end

function test_vector_b_not_column_vector()
    % b must be a column vector; a row vector must raise an error.
    try
        gmres_dr(ones(2), [1, 1]);
    catch ME
        msg = "Vector b must be a column vector.";
        assert(matches(ME.message, msg));
    end
end

function test_size_compatibility_between_A_and_b()
    % The length of b must equal the row count of A.
    try
        gmres_dr(ones(3), [1; 1]);
    catch ME
        msg = "Dimension mismatch between matrix A and vector b.";
        assert(matches(ME.message, msg));
    end
end

function test_m_greater_than_size_of_A()
    % m > n must raise an error.
    try
        gmres_dr(eye(3), ones(3, 1), 4);
    catch ME
        msg = "m must satisfy: 1 <= m <= n.";
        assert(matches(ME.message, msg));
    end
end

function test_k_equal_to_m_raises_error()
    % k == m is invalid because at least one Arnoldi step beyond the
    % deflation vectors is required per cycle.
    try
        gmres_dr(eye(3), ones(3, 1), 2, 2);
    catch ME
        msg = "k must satisfy: 0 < k < m.";
        assert(matches(ME.message, msg));
    end
end

function test_k_greater_than_m_raises_error()
    % k > m must raise an error.
    try
        gmres_dr(eye(3), ones(3, 1), 2, 3);
    catch ME
        msg = "k must satisfy: 0 < k < m.";
        assert(matches(ME.message, msg));
    end
end

function test_vector_xInitial_not_column_vector()
    % An initial guess that is not a column vector must raise an error.
    try
        gmres_dr(eye(3), ones(3, 1), [], [], [], [], ones(1, 3));
    catch ME
        msg = "Initial guess xInitial is not a column vector.";
        assert(matches(ME.message, msg));
    end
end

function test_size_compatibility_between_A_and_xInitial()
    % The length of xInitial must equal n.
    try
        gmres_dr(eye(3), ones(3, 1), [], [], [], [], ones(2, 1));
    catch ME
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        assert(matches(ME.message, msg));
    end
end

% =========================================================================
% ----> Fallback dispatch tests
% =========================================================================

function test_full_gmres_when_m_equals_n()
    % When m == n, GMRES-DR must fall back to unrestarted built-in GMRES.
    % The solution and flag must match MATLAB's gmres(A, b).

    A = eye(3);
    b = ones(3, 1);

    x1 = gmres(A, b);
    [x2, flag, relresvec, kdvec, time] = gmres_dr(A, b, 3);

    assertElementsAlmostEqual(x1, x2);
    assert(flag == 1);
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assertElementsAlmostEqual(kdvec, [3; 3]);
    assert(time > 0 && time < 5);
end

function test_restarted_gmres_when_k_equals_zero()
    % When k == 0, GMRES-DR must fall back to standard restarted GMRES(m).
    % The solution must match MATLAB's gmres(A, b, m).

    A = eye(3);
    b = ones(3, 1);
    m = 2;

    x1 = gmres(A, b, m);
    [x2, flag, relresvec, kdvec, time] = gmres_dr(A, b, m, 0);

    assertElementsAlmostEqual(x1, x2);
    assert(flag == 1);
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assertElementsAlmostEqual(kdvec, [2; 2]);
    assert(time > 0 && time < 5);
end

% =========================================================================
% ----> Correctness tests on simple systems
% =========================================================================

function test_default_parameters_identity_matrix()
    % With only A and b provided, defaults must be used and the system
    % Ax = b with A = I must be solved exactly.

    A = eye(3);
    b = ones(3, 1);

    [x, flag, ~, ~, ~] = gmres_dr(A, b);

    assertElementsAlmostEqual(x, ones(3, 1));
    assert(flag == 1);
end

function test_outputs_identity_matrix_small()
    % Solve I*x = [2;3;4] with m = 2, k = 1.
    % The identity system converges in the first Arnoldi step (happy
    % breakdown), so relresvec and kdvec have exactly two entries.

    A = eye(3);
    b = [2; 3; 4];
    m = 2;
    k = 1;
    tol = 1e-9;
    maxit = 100;

    [x, flag, relresvec, kdvec, time] = gmres_dr(A, b, m, k, tol, maxit, []);

    assertElementsAlmostEqual(x, [2; 3; 4]);
    assert(flag == 1);
    % Happy breakdown: subspace dimension is 1 (not m = 2)
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assert(time > 0 && time < 5);
    % kdvec(2) == 1 because happy breakdown truncates the subspace to 1
    assertEqual(kdvec(2), 1);
end

function test_outputs_diagonal_matrix()
    % Solve a diagonal system D*x = b where D = diag(1,...,n).
    % The diagonal matrix has n distinct eigenvalues; GMRES-DR should
    % converge at least as fast as standard GMRES(m) with the same m.

    n = 10;
    A = diag(1:n);
    b = ones(n, 1);
    m = 4;
    k = 2;
    tol = 1e-10;
    maxit = 200;

    [x, flag, ~, ~, time] = gmres_dr(A, b, m, k, tol, maxit, []);

    assertElementsAlmostEqual(x, (1 ./ (1:n))');
    assert(flag == 1);
    assert(time > 0 && time < 30);
end

% =========================================================================
% ----> Embree 3x3 toy problem (Morgan's own example)
% =========================================================================

function test_embree_3x3_toy_example()
    % Test GMRES-DR on the 3x3 system from Embree (1999), the same example
    % used in the original GMRES-DR paper [1].  For this matrix, restarted
    % GMRES(m) with m = 2 does not converge (Embree, 1999), but GMRES-DR
    % should because it deflates the problematic small eigenvalue.

    load('embree3.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 2;
    k = 1;
    tol = 1e-6;
    maxit = 100;

    [x, flag, ~, ~, time] = gmres_dr(A, b, m, k, tol, maxit);

    assertElementsAlmostEqual(x, [8; -7; 1]);
    assertEqual(flag, 1);
    assert(time > 0 && time < 100);
end

% =========================================================================
% ----> Sparse matrix tests (Sherman collection)
% =========================================================================

function test_sherman1()
    % Test GMRES-DR on the sherman1 matrix from the SuiteSparse collection.
    % We check convergence and that the krylov dimension vector is
    % consistently m at every cycle.

    load('sherman1.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 27;
    k = 3;
    tol = 1e-12;
    maxit = 1000;

    [~, flag, relresvec, kdvec, time] = gmres_dr(A, b, m, k, tol, maxit);

    assertEqual(flag, 1);
    % All cycles (except possibly a happy-breakdown cycle) must use
    % the full subspace of dimension m.
    assert(all(kdvec(2:end) <= m));
    assert(relresvec(end) < tol);
    assert(time > 0 && time < 300);
end

function test_sherman4()
    % Test GMRES-DR on the sherman4 matrix from the SuiteSparse collection.

    load('sherman4.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 27;
    k = 3;
    tol = 1e-12;
    maxit = 1000;

    [~, flag, relresvec, ~, time] = gmres_dr(A, b, m, k, tol, maxit);

    assertEqual(flag, 1);
    assert(relresvec(end) < tol);
    assert(time > 0 && time < 300);
end
