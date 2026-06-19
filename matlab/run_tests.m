function run_tests()
    %
    %   run tests with code coverage
    %
    %   USAGE::
    %
    %       run_tests()
    %
    %   [Modified from Remi Gau 2022 under GPL v3 license]
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

    tic;

    matlab_dir = fileparts(mfilename('fullpath'));
    cd(matlab_dir);

    fprintf('\nHome is %s\n', getenv('HOME'));

    folder_to_cover = fullfile(matlab_dir, 'src');

    test_folder = fullfile(matlab_dir, 'tests');

    addpath(genpath(folder_to_cover));
    addpath(fullfile(matlab_dir, '..', 'data'));

    if ispc
        success = moxunit_runtests(test_folder, '-verbose');

    else
        success = moxunit_runtests(test_folder, ...
                                   '-verbose', ...
                                   '-recursive', ...
                                   '-with_coverage', ...
                                   '-cover', ...
                                   folder_to_cover, ...
                                   '-cover_xml_file', ...
                                   'coverage.xml', ...
                                   '-cover_html_dir', ...
                                   fullfile(pwd, 'coverage_html') ...
                                  );
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
