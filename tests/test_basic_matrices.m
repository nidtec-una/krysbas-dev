function test_suite = test %#ok<*STOUT>
    %
    % Modified from:
    % https://github.com/Remi-Gau/template_matlab_analysis/blob/main/tests/test_my_fibonacci.m

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end

    initTestSuite;

end

function test_identity_matrix()
    % Test algorithms with identity matrix and rhs of ones.

    A = eye(3);
    b = ones(3, 1);
    x = pd_GMRES()

end