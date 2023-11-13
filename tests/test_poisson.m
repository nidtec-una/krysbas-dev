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
    
    % Source term and boundary values
    f = @(x) -6*x;
    g = @(x) x.^3;
    NODES = 9;
    INNERNODES = NODES-2;
    aStart = 0;
    aEnd = 1;
    h = ( aEnd - aStart ) / (NODES - 1);
    
    % Left-hand side (LHS) matrix 'A'
    % LHS is composed by the finite-difference coefficients for u''(x)
    A = 2*eye(INNERNODES) - diag(ones(INNERNODES-1,1), 1) - diag(ones(INNERNODES-1, 1), -1);
    
    % Right-hand side 'b' (RHS)
    % RHS is composed by the contributions from 
    % the boundary conditions g(x) and the source term f(x).
    b = h*h*f( linspace(aStart + h, aEnd - h, INNERNODES) )'; % b = h^2*f(x)
    % Add contribution from Dirichlet nodes to vector 'b'
    % b(1, 1) = g( astart );
    b(INNERNODES,1) = b(INNERNODES,1) + g( aEnd );
    
    % Exact solution 'uExact'
    uExact = g( linspace(aStart + h, aEnd - h, INNERNODES) )';
    
    % CALL ALGORITHMS: PD_GMRES, ADAPTIVE_GMRES, SWITCH_GMRES, etc.
    
    % PD_GMRES
    tol=1e-9; maxit=100; mInitial=3;
    uPD_GMRES= pd_gmres(A, b, mInitial, [], [], tol, maxit);

    % Assert whether the pd_gmres solution match the exact knonw solution
    assertElementsAlmostEqual(uPD_GMRES, uExact)

end