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


% ----> Test for sanity checks

function test_number_of_input_arguments() 
% Test if error is raised when passing incorrect number of inputs.  
    % Error should be raised since 1 parameter is given
    try
        pd_gmres(ones(2, 2));
    catch ME
        msg = "Too few input parameters. Expected at least A and b.";
        assert(ME.message == msg)
    end

    % Error should be raised since 10 parameters are given
    try 
        pd_gmres([], [], [], [], [], [], [], [], [], []);
    catch ME
        msg = "Too many input parameters.";
        assert(ME.message == msg)
    end     
end