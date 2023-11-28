function test_suite = test_pd_gmres %#ok<*STOUT>
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

    % Error should be raised since 10 parameters are given
    try
        pd_gmres([], [], [], [], [], [], [], [], [], []);
    catch ME
        msg = "Too many input parameters.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_matrix_A()
    % Test if error is raised when an empty matrix A is given.
    try
        pd_gmres([], 1);
    catch ME
        msg = 'Matrix A cannot be empty.';
        assert(matches(ME.message, msg));
    end
end

function test_non_square_matrix_A()
    % Test if error is raised when matrix A is not squared
    try
        pd_gmres([1; 1], [1; 1]);
    catch ME
        msg = "Matrix A must be square.";
        assert(matches(ME.message, msg));
    end
end

function test_empty_vector_b()
    % Test if error is raised when vector b is empty
    try
        pd_gmres(1, []);
    catch ME
        msg = "Vector b cannot be empty.";
    end
end

function test_vector_b_not_column_vector()
    % Test if error is raised when vector b is not a column vector
    try
        pd_gmres(ones(2), [1, 1]);
    catch ME
        msg = "Vector b must be a column vector.";
    end
end

function test_size_compatibility_between_A_and_b()
    % Test if error is raised when the dimensionality of A and b differs
    try
        pd_gmres(ones(3), [1; 1]);
    catch ME
        msg = "Dimension mismatch between matrix A and vector b.";
    end
end

function test_mInitial_not_given_empty_or_equal_to_n_use_unrestarted()
    % Test whether the unrestarted version of pd_gmres() is used if mInitial is
    % empty, not given, or equal to the dimension of the linear system n

    % Here, we actually solve a linear system and check whether
    % 'restarted' is false or not. This is not ideal, but it might be the
    % only way to construct the test

    A = eye(3);
    b = ones(3, 1);

    % Only provide A and b
    [~, ~, ~, ~, mvec, ~] = pd_gmres(A, b);
    assert(isnan(mvec));

    % Pass mInitial empty
    [~, ~, ~, ~, mvec, ~] = pd_gmres(A, b, []);
    assert(isnan(mvec));

    % Pass mInitial = 3
    [~, ~, ~, ~, mvec, ~] = pd_gmres(A, b, 3);
    assert(isnan(mvec));

end

function test_mInitial_given_use_restarted()
    % Test whether the restarted version of pd_gmres() is used if a valid
    % mInitial is provided

    A = eye(3);
    b = ones(3, 1);
    [~, ~, ~, ~, mvec, ~] = pd_gmres(A, b, 1);
    assert(~all(isnan(mvec)));

end

function test_mInitial_valid_range()
    % Test whether the initial restart parameter 'm' lies between (1, n).

    A = eye(4);
    b = ones(4, 1);

    % Loop over values that lie outside the range
    mInitialValues = [-50, 0, length(b) + 1, 100];
    for i = 1:length(mInitialValues)
        try
            pd_gmres(A, b, mInitialValues(i));
        catch ME
            msg = "mInitial must satisfy: 1 <= mInitial <= n.";
            assert(matches(ME.message, msg));
        end
    end

end

function test_warning_raised_if_mMinMax_given_when_unrestarted()
    % Test whether a warning is raised if mMinMax is given but restarted=false

    % Inputs that will generated the expected warning
    A = eye(3);
    b = ones(3, 1);
    mInitial = [];
    mMinMax = [1; 2];

    lastwarn('');  % Make sure to clear the last warning message
    warning('off');  % Avoid showing all warnings
    pd_gmres(A, b, mInitial, mMinMax);  % Call pd_gmres
    warning('on');  % Show all warnings again
    [warnMsg, ~] = lastwarn;  % retrieve warning message
    assert(matches(warnMsg, "mMinMax was given but will not be used."));

end

function test_mMinMax_valid_range()
    % Test the range of validity of the minimum and maximum values of m

    A = eye(3);
    b = ones(3, 1);
    mInitial = 1;

    mMinMaxValues = cell(3, 1);
    mMinMaxValues{1} = [0; 1];  % invalid mMin
    mMinMaxValues{2} = [1; 4];  % invalid mMax
    mMinMaxValues{3} = [1; 1];  % non-monotonically increasing

    msg = "mMinMax must satisfy: 1 <= mMinMax(1) < mMinMax(2) <= n.";
    for ii = 1:length(mMinMaxValues)
        try
            pd_gmres(A, b, mInitial, mMinMaxValues{ii});
        catch ME
            assert(matches(ME.message, msg));
        end
    end

end

function test_mMinMax_valid_range_wrt_mInitial()
    % Test the validity of mMinMax with respect to mInitial

    A = eye(3);
    b = ones(3, 1);
    mInitial = 1;
    mMinMax = [2; 3];

    try
        pd_gmres(A, b, mInitial, mMinMax);
    catch ME
        msg = 'mMinMax must satisfy: mMinMax(1) <= mInitial <= mMinMax(2).';
        assert(matches(ME.message, msg));
    end

end

function test_warning_raised_if_mStep_given_when_unrestarted()
    % Test if a warning is raised when mStep is given but restarted=false

    % Inputs that will generated the expected warning
    A = eye(3);
    b = ones(3, 1);
    mInitial = [];
    mMinMax = [];
    mStep = 2;

    lastwarn('');  % Make sure to clear the last warning message
    warning('off');  % Avoid showing all warnings
    pd_gmres(A, b, mInitial, mMinMax, mStep);  % Call pd_gmres
    warning('on');  % Show all warnings again
    [warnMsg, ~] = lastwarn;  % retrieve warning message
    assert(matches(warnMsg, "mStep was given but will not be used."));

end

function test_mStep_valid_range()
    % Test valid range of parameter mStep

    A = eye(3);
    b = ones(3, 1);
    mInitial = 1;
    mMinMax = [1; 3];

    mStepValues = [0, 3, 10];  % these values lie outside the valid range

    for ii = 1:length(mStepValues)
        try
            pd_gmres(A, b, mInitial, mMinMax, mStepValues(ii));
        catch ME
            msg = "mStep must satisfy: 0 < mStep < n.";
            assert(matches(ME.message, msg));
        end
    end

end

function test_vector_xInitial_not_column_vector()
    % Test if error is raised when xInitial is not a column vector

    A = eye(3);
    b = ones(3, 1);
    x0 = ones(1, 3);

    try
        pd_gmres(A, b, [], [], [], [], [], x0);
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
        pd_gmres(A, b, [], [], [], [], [], x0);
    catch ME
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        assert(matches(ME.message, msg));
    end
end

function test_outputs_unrestarted_identity_matrix()
    % Test whether the correct outputs are returned using an identity
    % matrix when the unrestarted algorithm is called

    % Setup a trivial linear system
    A = eye(3);
    b = ones(3, 1);
    
    % Call function                                             
    [x, flag, relres, resvec, mvec, time] = pd_gmres(A, b);
   
    % Compare with expected outputs
    assertEqual(x, ones(3, 1));
    assert(flag == 1)
    assertEqual(relres, 0)
    assertElementsAlmostEqual(resvec, [1.7320508075688770; 0])
    assert(isnan(mvec))
    assert(time > 0 && time < 5) 
end

function test_outputs_restarted_identity_matrix()
    % Test whether the correct outputs are returned using an identity
    % matrix when the restarted algorithm is called

    A = eye(3);
    b = [2; 3; 4];
    mInitial = 1;
    mMinMax = [1; 3];
    mStep = 1;
    tol = 1e-9;
    maxit = 10;
    xInitial = zeros(3, 1);
    alphaPD = [-3; 5];
       
    % Call function
    [x, flag, relres, resvec, mvec, time] = ...
        pd_gmres(A, b, mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD);

    % Compare with expected outputs
    assertElementsAlmostEqual(x, [2; 3; 4]);
    assert(flag == 1)
    assert(relres <= tol)
    assertElementsAlmostEqual(resvec, [1; 0])
    assertEqual(mvec, [1; 1])
    assert(time > 0 && time < 5) 

end