function run_tests()
    %
    % run tests with code coverage
    %
    % USAGE::
    %
    %   run_tests()
    %
    % Modified from Remi Gau 2022 under GPL v3 license.

    tic;

    cd(fileparts(mfilename('fullpath')));

    fprintf('\nHome is %s\n', getenv('HOME'));

    folder_to_cover = fullfile(pwd, 'src');

    test_folder = fullfile(pwd, 'tests');

    addpath(genpath(folder_to_cover));

    if ispc
        success = moxunit_runtests(test_folder, '-verbose');

    else
        success = moxunit_runtests(test_folder, ...
                                   '-verbose', '-recursive', '-with_coverage', ...
                                   '-cover', folder_to_cover, ...
                                   '-cover_xml_file', 'coverage.xml', ...
                                   '-cover_html_dir', fullfile(pwd, 'coverage_html'));
    end

    fileID = fopen('test_report.log', 'w');
    if success
        fprintf(fileID, '0');
    else
        fprintf(fileID, '1');
    end
    fclose(fileID);

    toc;

end