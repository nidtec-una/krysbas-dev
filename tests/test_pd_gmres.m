function test_suite = test_pd_gmres %#ok<*STOUT>
    %
    % Modified from:
    % https://github.com/Remi-Gau/template_matlab_analysis/blob/main/tests/test_my_fibonacci.m

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end

    initTestSuite;

end


% ----> Test sanity checks and default values

function test_number_of_input_arguments() 
% Test if error is raised when passing incorrect number of inputs.  
    % Error should be raised since 1 parameter is given
    try
        pd_gmres(ones(2, 2));
    catch ME
        msg = "Too few input parameters. Expected at least A and b.";
        assert(matches(ME.message, msg))
    end

    % Error should be raised since 10 parameters are given
    try 
        pd_gmres([], [], [], [], [], [], [], [], [], []);
    catch ME
        msg = "Too many input parameters.";
        assert(matches(ME.message, msg))
    end     
end

function test_empty_matrix_A()
% Test if error is raised when an empty matrix A is given.
    try
        pd_gmres([], 1);
    catch ME
        msg = 'Matrix A cannot be empty.';
        assert(matches(ME.message, msg))
    end
end

function test_non_square_matrix_A()
% Test if error is raised when matrix A is not squared
    try
        pd_gmres([1; 1], [1; 1]);
    catch ME
        msg = "Matrix A must be square.";
        assert(matches(ME.message, msg))
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

function test_mInitial_not_given_empty_or_equal_to_n()
% Test whether the unrestarted version of pd_gmres() is used if mInitial is
% empty, not given, or equal to the dimension of the linear system n

    % Here, we actually solve a linear system and check whether
    % 'restarted' is false or not. This is not ideal, but it might be the
    % only way to construct the test

    A = eye(3);
    b = ones(3, 1);
    
    % Only provide A and b
    [~, ~, ~, ~, ~, restarted, ~] = pd_gmres(A, b);
    assert(restarted == false);

    % Pass mInitial empty
    [~, ~, ~, ~, ~, restarted, ~] = pd_gmres(A, b, []);
    assert(restarted == false);

    % Pass mInitial = 3
    [~, ~, ~, ~, ~, restarted, ~] = pd_gmres(A, b, 3);
    assert(restarted == false);

end

function test_mInitial_valid_range()
% Test whether the initial restart parameter 'm' lies between (1, n).

    A = eye(4);
    b = ones(4, 1);

    % Loop over values that lie outside the range
    mInitialValues = [-50, 0, length(b) + 1, 100];
    for i=1:length(mInitialValues)
        try
            pd_gmres(A, b, mInitialValues(i))
        catch ME
            msg = "mInitial must satisfy: 1 <= mInitial <= n.";
            assert(matches(ME.message, msg))
        end
    end
            
end

