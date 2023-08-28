function test_suite = test_pd_GMRES %#ok<*STOUT>
    %
    % Modified from:
    % https://github.com/Remi-Gau/template_matlab_analysis/blob/main/tests/test_my_fibonacci.m

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end

    initTestSuite;

end

function test_pd_GMRES_0()

    sol7=[]; %PD-GMRES(m)

    % Load A and b from the sherman3 matrix
    rootFolder = fileparts(pwd); % go to the root folder
    dataFolder = fullfile(rootFolder, 'data'); % enter the data folder
    load 'sherman3.mat' Problem  %

    % Matrix A and right-hand-side b
    A=Problem.A;
    b=Problem.b;

    % Parameters
    alpha=-3; %2
    delta=5; %0.8
    opts_tol=1e-9;
    itermax=1000;
    p = 1;

    %%       %PD-GMRES(m)
    % The original idea was to compute the average execution time,
    % we may discuss if this is still necessary 
    for i=1:p
        color_pd_gmres='b';
        mPD=30;
        %alpha=2;
        [time, logres_pd_gmres]=pd_GMRES(A,b, mPD, alpha, delta,itermax);
        sol7(size(sol7,1)+1,:)= [time];
    end

end