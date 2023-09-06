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

% function test_pd_GMRES_0()  
% 
%     sol7=[]; %PD-GMRES(m)
% 
%     % Load A and b from the sherman3 matrix
%     rootFolder = fileparts(pwd); % go to the root folder
%     dataFolder = fullfile(rootFolder, 'data'); % enter the data folder
%     exampleProblem = fullfile(dataFolder, 'sherman3.mat');
%     load(exampleProblem); %   Problem  %
% 
%     % Matrix A and right-hand-side b
%     A=Problem.A;
%     b=Problem.b;
% 
%     % Parameters
%     alpha=-3; %2
%     delta=5; %0.8
%     opts_tol=1e-9;
%     itermax=1000;
%     p = 1;
% 
%     %%       %PD-GMRES(m)
%     % The original idea was to compute the average execution time,
%     % we may discuss if this is still necessary 
%     for i=1:p
%         color_pd_gmres='b';
%         mPD=30;
%         %alpha=2;
%         rootFolder = fileparts(pwd); % go to the root folder
%         srcFolder = fullfile(rootFolder, 'src'); % enter the data folder
%         cd(srcFolder)
%         [time, logres_pd_gmres, xx]=pd_gmres(A,b, mPD, alpha, delta,itermax);
%         sol7(size(sol7,1)+1,:)= [time];
%     end
% 
% end

function test_pd_GMRES_1()
    % Solve a 1D-Poisson equation at [aStart, aEnd] = [0, 1]
    % u''(x) = f(x), u(0) = g(0) = 0, u(1) = g(1) = 1;
    
    % f(x) = 6*x; g(x) = x^3
    f = @(x) 6*x;
    g = @(x) x.^3;
    
    NODES = 41;
    aStart = 0;
    aEnd = 1;
    h = ( aEnd - aStart ) / (NODES - 1);
    b = h*h*f( linspace(aStart, aEnd, NODES) )';
    b(1,1) = g( aStart );
    b(NODES,1) = g( aEnd );
    
    A = 2*eye(NODES) - diag(ones(NODES-1,1), 1) - diag(ones(NODES-1, 1), -1);
    A(1, 1:2) = [1 0];
    A(NODES, NODES-1:NODES) = [0 1];
    
    u1 = A\b;
    u2 = g( linspace(aStart, aEnd, NODES) )';
    
    rootFolder = fileparts(pwd); % go to the root folder
    srcFolder = fullfile(rootFolder, 'src'); % enter the data folder
    cd(srcFolder)
    %     [time, logres_pd_gmres]=pd_gmres(A,b, mPD, alpha, delta,itermax); 
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
        rootFolder = fileparts(pwd); % go to the root folder
        srcFolder = fullfile(rootFolder, 'src'); % enter the data folder
        cd(srcFolder)
        [time, logres_pd_gmres, u3]=pd_gmres(A,b, mPD, alpha, delta,itermax);
        %sol7(size(sol7,1)+1,:)= [time];
    end
    
    % We probably need to add here a moxunit test assertElementsAlmostEqual(u1,u3)
    u4 = u3 - u1;
    disp(max(u4))

end
