function test_suite = test_slgmres_e %#ok<*STOUT>
    %
    %   Test suite for the SLGMRES-E solver.
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

    try
        slgmres_e(ones(2, 2));
    catch ME
        msg = "Too few input parameters. Expected at least A and b.";
        assert(matches(ME.message, msg));
    end

    try
        slgmres_e([], [], [], [], [], [], [], [], [], [], []);
    catch ME
        msg = "Too many input parameters.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_matrix_A()
    try
        slgmres_e([], ones(2, 1));
    catch ME
        msg = "Matrix A cannot be empty.";
        assert(matches(ME.message, msg));
    end
end

function test_non_square_matrix_A()
    try
        slgmres_e([1; 1], [1; 1]);
    catch ME
        msg = "Matrix A must be square.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_vector_b()
    try
        slgmres_e(eye(2), []);
    catch ME
        msg = "Vector b cannot be empty.";
        assert(matches(ME.message, msg));
    end
end

function test_vector_b_not_column_vector()
    try
        slgmres_e(ones(2), [1, 1]);
    catch ME
        msg = "Vector b must be a column vector.";
        assert(matches(ME.message, msg));
    end
end

function test_size_compatibility_between_A_and_b()
    try
        slgmres_e(ones(3), [1; 1]);
    catch ME
        msg = "Dimension mismatch between matrix A and vector b.";
        assert(matches(ME.message, msg));
    end
end

function test_m_greater_than_size_of_A()
    try
        slgmres_e(eye(3), ones(3, 1), 4);
    catch ME
        msg = "m must satisfy: 1 <= m <= n.";
        assert(matches(ME.message, msg));
    end
end

function test_l_not_positive_raises_error()
    try
        slgmres_e(eye(3), ones(3, 1), 2, 0);
    catch ME
        msg = "l must satisfy: l > 0.";
        assert(matches(ME.message, msg));
    end
end

function test_d_not_positive_raises_error()
    try
        slgmres_e(eye(3), ones(3, 1), 2, 1, 0);
    catch ME
        msg = "d must satisfy: d > 0.";
        assert(matches(ME.message, msg));
    end
end

function test_epsilon0_out_of_range_raises_error()
    try
        slgmres_e(eye(3), ones(3, 1), 2, 1, 1, 1.5);
    catch ME
        msg = "epsilon0 must satisfy: 0 < epsilon0 < 1.";
        assert(matches(ME.message, msg));
    end
end

function test_vector_xInitial_not_column_vector()
    try
        slgmres_e(eye(3), ones(3, 1), [], [], [], [], [], [], ones(1, 3));
    catch ME
        msg = "Initial guess xInitial is not a column vector.";
        assert(matches(ME.message, msg));
    end
end

function test_size_compatibility_between_A_and_xInitial()
    try
        slgmres_e(eye(3), ones(3, 1), [], [], [], [], [], [], ones(2, 1));
    catch ME
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        assert(matches(ME.message, msg));
    end
end

% =========================================================================
% ----> Fallback dispatch tests
% =========================================================================

function test_full_gmres_when_m_equals_n()
    % When m == n, SLGMRES-E must fall back to unrestarted built-in GMRES.

    A = eye(3);
    b = ones(3, 1);

    x1 = gmres(A, b);
    [x2, flag, relresvec, kdvec, time] = slgmres_e(A, b, 3);

    assertElementsAlmostEqual(x1, x2);
    assert(flag == 1);
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assertElementsAlmostEqual(kdvec, [3; 3]);
    assert(time > 0 && time < 5);
end

% =========================================================================
% ----> Correctness tests on simple systems
% =========================================================================

function test_default_parameters_identity_matrix()
    A = eye(3);
    b = ones(3, 1);

    [x, flag, ~, ~, ~] = slgmres_e(A, b);

    assertElementsAlmostEqual(x, ones(3, 1));
    assert(flag == 1);
end

function test_outputs_identity_matrix_small()
    % Solve I*x = [2;3;4] with m = 2, l = 1, d = 1. The identity system
    % converges on the very first (plain) cycle, before any switching
    % decision is exercised. A*v1 is parallel to v1 for A = I, so the
    % Arnoldi loop hits a happy breakdown after a single step regardless
    % of m -- kdvec(2) == 1, matching the same case documented in
    % test_gmres_dr.m.

    A = eye(3);
    b = [2; 3; 4];
    m = 2;
    l = 1;
    d = 1;
    tol = 1e-9;
    maxit = 100;

    [x, flag, relresvec, kdvec, time] = ...
        slgmres_e(A, b, m, l, d, [], tol, maxit);

    assertElementsAlmostEqual(x, [2; 3; 4]);
    assert(flag == 1);
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assertEqual(kdvec(2), 1);
    assert(time > 0 && time < 5);
end

% =========================================================================
% ----> Embree 3x3 toy problem
% =========================================================================

function test_embree_3x3_toy_example()
    % Test SLGMRES-E on the 3x3 system from Embree (1999). m = 2 leaves no
    % room for a distinct "fresh" subspace once either augmentation type
    % is added (m + l or m + d exceeds n only mildly), but the switching
    % logic itself is still exercised meaningfully across cycles.

    load('embree3.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 2;
    l = 1;
    d = 1;
    tol = 1e-6;
    maxit = 100;

    [x, flag, ~, ~, time] = slgmres_e(A, b, m, l, d, [], tol, maxit);

    assertElementsAlmostEqual(x, [8; -7; 1], 'relative', 1e-4);
    assertEqual(flag, 1);
    assert(time > 0 && time < 100);
end

% =========================================================================
% ----> Sparse matrix tests
% =========================================================================

function test_sherman1()
    load('sherman1.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 27;
    l = 3;
    d = 3;
    tol = 1e-12;
    maxit = 1000;

    [~, flag, relresvec, ~, time] = slgmres_e(A, b, m, l, d, [], tol, maxit);

    assertEqual(flag, 1);
    assert(relresvec(end) < tol);
    assert(time > 0 && time < 300);
end

function test_sherman4()
    load('sherman4.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 27;
    l = 3;
    d = 3;
    tol = 1e-12;
    maxit = 1000;

    [~, flag, relresvec, ~, time] = slgmres_e(A, b, m, l, d, [], tol, maxit);

    assertEqual(flag, 1);
    assert(relresvec(end) < tol);
    assert(time > 0 && time < 300);
end

function test_sherman5_matches_cabral_reference()
    % Regression test against Cabral, Schaerer & Bhaya (2020)'s own
    % reference implementation (Adaptive_lgmres_e_switch.m), run with the
    % exact parameters from their master_algoritmos.m driver. Verified by
    % direct numerical comparison: the reference converges in 363 cycles
    % to a final relative residual of 9.294759e-10; this port must match
    % essentially exactly (both use the same switching decision at every
    % single cycle -- verified cycle-by-cycle during development).

    load('sherman5.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    m = 28;
    l = 2;
    d = 2;
    epsilon0 = 0.01;
    tol = 1e-9;
    maxit = 1000;

    [~, flag, relresvec, ~, time] = ...
        slgmres_e(A, b, m, l, d, epsilon0, tol, maxit);

    assertEqual(flag, 1);
    assertEqual(length(relresvec), 364);
    assertElementsAlmostEqual(relresvec(end), 9.294759e-10, 'relative', 1e-2);
    assert(time > 0 && time < 300);
end
