function test_suite = test_lgmres %#ok<*STOUT>
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

% ----> Test for sanity checks

function test_number_of_input_arguments()
    % Test if error is raised when passing incorrect number of inputs.
    % Error should be raised since 1 parameter is given
    try
        pd_gmres(ones(2, 2));
    catch ME
        msg = "Too few input parameters. Expected at least A and b.";
        assert(matches(ME.message, msg));
    end

    % Error should be raised since 8 parameters are given
    try
        lgmres([], [], [], [], [], [], [], []);
    catch ME
        msg = "Too many input parameters.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_matrix_A()
    % Test if error is raised when an empty matrix A is given.
    try
        lgmres([], 1);
    catch ME
        msg = 'Matrix A cannot be empty.';
        assert(matches(ME.message, msg));
    end
end

function test_non_square_matrix_A()
    % Test if error is raised when matrix A is not squared
    try
        lgmres([1; 1], [1; 1]);
    catch ME
        msg = "Matrix A must be square.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_vector_b()
    % Test if error is raised when vector b is empty
    try
        lgmres(1, []);
    catch ME
        msg = "Vector b cannot be empty.";
    end
end

function test_vector_b_not_column_vector()
    % Test if error is raised when vector b is not a column vector
    try
        lgmres(ones(2), [1, 1]);
    catch ME
        msg = "Vector b must be a column vector.";
    end
end

function test_size_compatibility_between_A_and_b()
    % Test if error is raised when the dimensionality of A and b differs
    try
        lgmres(ones(3), [1; 1]);
    catch ME
        msg = "Dimension mismatch between matrix A and vector b.";
    end
end

% Here we write tests for m and k

function test_vector_xInitial_not_column_vector()
    % Test if error is raised when xInitial is not a column vector

    A = eye(3);
    b = ones(3, 1);
    x0 = ones(1, 3);

    try
        lgmres(A, b, [], [], [], [], x0);
    catch ME
        msg = "Initial guess xInitial is not a column vector.";
        assert(matches(ME.message, msg));
    end
end

function test_size_compatibility_between_A_and_xInitial()
    % Test size compatibility between A and xInitial

    A = eye(3);
    b = ones(3, 1);
    x0 = ones(2, 1);

    try
        lgmres(A, b, [], [], [], [], x0);
    catch ME
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        assert(matches(ME.message, msg));
    end
end

% Check the name of this function
function test_outputs_unrestarted_identity_matrix() % Linear system # 1
    % Test whether the correct outputs are returned using an identity
    % matrix when the unrestarted algorithm is called

    % Setup a trivial linear system
    n = 100;
    A = eye(n);
    b = ones(n, 1);
    xInitial = zeros(n, 1);
    
    % Call LGMRES
    [x, flag, relresvec, time] = lgmres(A, b, 27, 3, 1e-6, 100, xInitial);
    
    % Compare with expected outputs
    assertElementsAlmostEqual(x, ones(n, 1));
    assert(flag == 1);
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assert(time > 0 && time < 5);
end

% Check the name of this function
function test_outputs_restarted_identity_matrix() % Linear system # 2
    % Test whether the correct outputs are returned using an identity
    % matrix when the restarted algorithm is called

    % Setup a trivial linear system
    n = 3;
    A = eye(n);
    b = [2; 3; 4];
    xInitial = zeros(n, 1);
    
    % Call LGMRES
    [x, flag, relresvec, time] = lgmres(A, b, 2, 1, 1e-6, 100, xInitial);
    
    % Compare with expected outputs
    assertElementsAlmostEqual(x, [2; 3; 4]);
    assert(flag == 1);
    assertElementsAlmostEqual(relresvec, [1; 0]);
    assert(time > 0 && time < 5);

end

function test_embree_three_by_three_toy_example()
    % Test Embree's 3x3 linear system from https://www.jstor.org/stable/25054403
    % with PD-GMRES. We construct two checks, one with mInitial = 1, and
    % one with mInitial = 2. Note that the gmres(m) algorithm with m = 2,
    % does not converge (this was proven by Embree in his paper). This, however,
    % does not occur with the PD-GMRES, which adaptively changes the value
    % of m to avoid stagnation.

    % Load A and b
    load('data/embree3.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    % Setup PD-GMRES
    mInitials = [1; 2];
    tol = 1e-9;
    maxit = 20;

    % Loop over the values of mInitials, i.e., {1, 2}.
    for i = 1:length(mInitials)
        % Call PD-GMRES
        [x, flag, relresvec, mvec, ~] = ...
            pd_gmres(A, b, mInitials(i), [], [], tol, maxit, [], []);

        % Compare with expected outputs
        assertElementsAlmostEqual(x, [8; -7; 1]);
        assertEqual(flag, 1);
        if mInitials(i) == 1
            relresvecExpected = [1.000000000000000
                                 0.925820099772551
                                 0.654653670707977
                                 0];
            mvecExpected = [1; 1; 1; 2];
        else
            relresvecExpected = [1.000000000000000
                                 0.462910049886276
                                 0.377189160453301
                                 0.000000000000001];
            mvecExpected = [2; 2; 2; 3];
        end
        assertElementsAlmostEqual(relresvec, relresvecExpected);
        assertEqual(mvec, mvecExpected);
    end

end

function test_sherman_one()
    % Test pd_gmres with sherman1 matrix

    % Load A and b
    load('data/sherman1.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    % Setup PD-GMRES
    mInitial = 30;
    tol = 1e-9;
    mStep = 3;
    maxit = 1000;

    % Call PD-GMRES
    [~, flag, ~, mvec, ~] = ...
            pd_gmres(A, b, mInitial, [], mStep, tol, maxit, [], []);

    % We check if it has converged and the total sum of outer iterations
    assertEqual(flag, 1);
    assertEqual(sum(mvec), 955);

end

function test_sherman_four()
    % Test pd_gmres with sherman4 matrix

    % Load A and b
    load('data/sherman4.mat', 'Problem');
    A = Problem.A;
    b = Problem.b;

    % Setup PD-GMRES
    mInitial = 30;
    tol = 1e-9;
    mStep = 3;
    maxit = 1000;

    % Call PD-GMRES
    [~, flag, ~, mvec, ~] = ...
            pd_gmres(A, b, mInitial, [], mStep, tol, maxit, [], []);

    % We check if it has converged and the total sum of outer iterations
    assertEqual(flag, 1);
    assertEqual(sum(mvec), 440);

end
