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

function test_pd_GMRES_01_poisson() 
    % Change name (change test?) by boundary conditions 
    % (Dirichlet, Neumann, Robin, mixed).
    % Description:
    % ============
    % Solve a 1D-Poisson equation at [aStart, aEnd] = [0, 1]
    % - u''(x) = f(x)
    % subject to
    % u(0) = g(0) = 0, u(1) = g(1) = 1;
    % where
    % f(x) = 6*x; g(x) = x^3
    % Discretization results in a linear system Au = b
    
    % Discretization grid
    f = @(x) -6*x;
    g = @(x) x.^3;
    NODES = 9;
    INNERNODES = NODES-2;
    aStart = 0;
    aEnd = 1;
    h = ( aEnd - aStart ) / (NODES - 1);
    
    % Left-hand side matrix 'A'
    A = 2*eye(INNERNODES) - diag(ones(INNERNODES-1,1), 1) - diag(ones(INNERNODES-1, 1), -1);
    
    % Right-hand side 'b'
    b = h*h*f( linspace(aStart + h, aEnd - h, INNERNODES) )'; % b = h^2*f(x)
    % Dirichlet conditions are here
    b(INNERNODES,1) = b(INNERNODES,1) + g( aEnd );
    
    % Exact solution 'uExact'
    uExact = g( linspace(aStart + h, aEnd - h, INNERNODES) )';
    
    % CALL ALGORITHMS: PD_GMRES, ADAPTIVE_GMRES, SWITCH_GMRES, etc.
    % PD_GMRES(m)
    % The original idea was to compute the average execution time,
    % we may discuss if this is still necessary 

    % Parameters for PD_GMRES
    alpha=-3;
    delta=5;
    opts_tol=1e-9;
    itermax=1000;
    p = 1;
    for i=1:p
        color_pd_gmres='b';
        mPD=3;
        uBackslash  = A\b;
        [~, uPD_GMRES]=pd_GMRES(A,b, mPD, alpha, delta,itermax,opts_tol);
        uExact = g( linspace(aStart + h, aEnd - h, INNERNODES) )';
    end

    assertElementsAlmostEqual(uPD_GMRES, uExact)

end