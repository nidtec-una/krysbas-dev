function test_suite = test_my_fibonacci %#ok<*STOUT>
    %
    % (C) Copyright 2022 Remi Gau
    % Adapted from Remi Gau 22

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end

    initTestSuite;

end

function test_my_summa()

    % THEN
    % assertEqual(results, [0, 1, 1, 2, 3, 5, 8, 13]);
    assertEqual(2, 1 + 1)

end

%function test_my_fibonacci_default()
%
%    % WHEN
%    results = my_fibonacci();
%
%    % THEN
%    assertEqual(results, [0, 1, 1, 2, 3, 5, 8]);
%
%end
%
%function test_my_fibonacci_error()
%
%    if ~is_octave()
%        assertExceptionThrown(@()my_fibonacci(-1),
%'MATLAB:InputParser:ArgumentFailedValidation');
%    end
%
%end