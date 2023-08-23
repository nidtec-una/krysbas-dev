% run tests with code coverage via the run_tests scripts in the root folder.
% Modified under GPL V3 License from:
% https://github.com/Remi-Gau/template_matlab_analysis/blob/main/.github/workflows/run_tests_ci.m

root_dir = getenv('GITHUB_WORKSPACE');

% MOxUnit and MOcov need to be in the matlab path
addpath(fullfile(root_dir, 'MOcov', 'MOcov'));  % check if this works
cd(fullfile(root_dir, 'MOxUnit', 'MOxUnit'));
run moxunit_set_path();

cd(root_dir);
run run_tests();