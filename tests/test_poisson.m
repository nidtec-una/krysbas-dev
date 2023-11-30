function test_suite = test_poisson %#ok<*STOUT>
    %
    %   Modified from:
    %   https://github.com/Remi-Gau/template_matlab_analysis
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

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end

    initTestSuite;

end

function test_poisson_1d_dir_bc()
    % Test solvers for a 1D-FDM linear system with Dirichlet boundaries.
    %
    % Description:
    % ============
    %   Solve a 1D-Poisson equation in the domain [aStart, aEnd] = [0, 1]
    %
    %                         - u''(x) = f(x)
    %
    %   subject to
    %
    %   	        u(0) = g(0) = 0, u(1) = g(1) = 1; 
    % 
    %   where 
    %   
    %           f(x) = 6*x          and         g(x) = x^3.
    %
    %   Discretization results in a linear system Au = b

    % Source term and boundary values
    f = @(x) -6 * x;
    g = @(x) x.^3;
    NODES = 9;
    INNERNODES = NODES - 2;
    aStart = 0;
    aEnd = 1;
    h = (aEnd - aStart) / (NODES - 1);

    % Left-hand side (LHS) matrix 'A'
    % LHS is composed by the finite-difference coefficients for u''(x)
    A = 2 * eye(INNERNODES) - diag(ones(INNERNODES - 1, 1), 1) - ...
        diag(ones(INNERNODES - 1, 1), -1);

    % Right-hand side 'b' (RHS)
    % RHS is composed by the contributions from
    % the boundary conditions g(x) and the source term f(x).
    b = h * h * f(linspace(aStart + h, aEnd - h, INNERNODES))'; % b = h^2*f(x)
    % Add contribution from Dirichlet nodes to vector 'b'
    b(INNERNODES, 1) = b(INNERNODES, 1) + g(aEnd);

    % Exact solution 'uExact'
    uExact = g(linspace(aStart + h, aEnd - h, INNERNODES))';

    % CALL ALGORITHMS: PD_GMRES, ADAPTIVE_GMRES, SWITCH_GMRES, etc.

    % PD_GMRES
    tol = 1e-9;
    maxit = 100;
    mInitial = 3;
    uPD_GMRES = pd_gmres(A, b, mInitial, [], [], tol, maxit);

    % Assert whether the pd_gmres solution match the exact knonw solution
    assertElementsAlmostEqual(uPD_GMRES, uExact);

end

function test_poisson_1d_dir_neu_bc()
    % Test solvers for a 1D-FDM linear system with 
    % Dirichlet-Neumann boundaries.
    %
    % Description:
    % ============
    %   Solve a 1D-Poisson equation in the domain [aStart, aEnd] = [0, 1]
    %
    %                         - u''(x) = f(x)
    %
    %   subject to
    %
    %                   u(0) = g(0) = 1, u'(1) = 4;
    %
    %   where
    %
    %            f(x) = - 2       and      g(x) = ( x + 1 )^2.
    %
    % \Discretization results in a linear system Au = b

    % Source term and boundary values
    f = - 2;
    g = @(x) ( x + 1 ).^2;
    NODES = 5;
    INNERNODES = NODES - 2;
    aStart = 0;
    aEnd = 1;
    h = (aEnd - aStart) / (NODES - 1);
    fx = -2; % for all the domain
    DgN= 4; % Neumann boundary constraint.

    % Left-hand side (LHS) matrix 'A'
    % LHS is composed by the finite-difference coefficients for u''(x)
    A = 2 * eye(INNERNODES+1) - diag(ones(INNERNODES+1 - 1, 1), 1) - ...
        diag(ones(INNERNODES+1 - 1, 1), -1);
    % Contribution due to Neumann boundary condition
    A(INNERNODES+1, INNERNODES+1) = 1;

    % Right-hand side 'b' (RHS)
    % RHS is composed by the contributions from
    % the boundary conditions g(x) and the source term f(x).
    % In this example, f(x) = f = -2 for all x.
    b = h * h * f * ones(INNERNODES+1, 1); % b = h^2*f(x)
    % Contribution from Dirichlet-Neumann nodes to vector 'b'
    % Dirichlet boundary condition at x = 0
    b(1, 1) = b(1, 1) + g(aStart); 
    % Neumann boundary condition at x = 1
    b(INNERNODES+1, 1) = 0.5 * h * h * f + h * DgN; 

    % Exact solution 'uExact'
    uExact = g(linspace(aStart + h, aEnd, INNERNODES+1))';

    % CALL ALGORITHMS: PD_GMRES, ADAPTIVE_GMRES, SWITCH_GMRES, etc.

    % PD_GMRES
    tol = 1e-9;
    maxit = 100;
    mInitial = 3;
    uPD_GMRES = pd_gmres(A, b, mInitial, [], [], tol, maxit);

    % Assert whether the pd_gmres solution match the exact knonw solution
    assertElementsAlmostEqual(uPD_GMRES, uExact);

end

function test_poisson_1d_robin_bc()
    % Test solvers for a 1D-FDM linear system with 
    % Robin-Robin boundaries.
    %
    % Description:
    % ============
    %   Solve a 1D-Poisson equation in the domain [aStart, aEnd] = [0, 1]
    %
    %                         u''(x) = f(x)
    %
    %   subject to Robin-Robin boundary conditions
    %           
    %       	        alpha * u(0) + u'(0) = 3, 
    %                   beta * u(1) + u'(1) = 0,
    %
    %   where   
    % 
    %           f(x) = 2        and        g(x) = ( x + 1 )^2.
    %
    %   Discretization results in a linear system Au = b

    % Source term and boundary values
    f = 2;
    g = @(x) ( x + 1 ).^2;
    NODES = 5;
    aStart = 0;
    aEnd = 1;
    alpha = 1; 
    beta = -1;
    cStart = 3;
    cEnd = 0;
    h = (aEnd - aStart) / (NODES - 1);
    fx = -2; % for all the domain
    DgN= 4; % Neumann boundary constraint.

    % Left-hand side (LHS) matrix 'A'
    % LHS is composed by the finite-difference coefficients for u''(x)
    A = full(gallery('tridiag', NODES, 1, -2, 1));
    % Contribution due to Robin boundary condition
    A(1, 1) = h * alpha - 1;
    A(NODES, NODES) = - h * beta - 1;

    % Right-hand side 'b' (RHS)
    % RHS is composed by the contributions from
    % the boundary conditions g(x) and the source term f(x).
    % In this example, f(x) = f = 2 for all x.
    b = h * h * f * ones(NODES, 1); % b = h^2*f(x)
    % Contribution from Robin-Robin nodes to vector 'b'
    % Robin boundary condition at x = 0
    b(1, 1) = 0.5 * h * h * f + h * cStart; 
    % Robin boundary condition at x = 1
    b(NODES, 1) = 0.5 * h * h * f - h * cEnd; 

    % Exact solution 'uExact'
    uExact = g(linspace(aStart, aEnd, NODES))';

    % CALL ALGORITHMS: PD_GMRES, ADAPTIVE_GMRES, SWITCH_GMRES, etc.

    % PD_GMRES
    tol = 1e-9;
    maxit = 100;
    mInitial = 3;
    uPD_GMRES = pd_gmres(A, b, mInitial, [], [], tol, maxit);

    % Assert whether the pd_gmres solution match the exact knonw solution
    assertElementsAlmostEqual(uPD_GMRES, uExact);

end

function test_dir_2d_dir_bc()
    % Test solvers for a 2D-FDM linear system with 
    % Dirichlet boundaries.
    %
    % Description:
    % ============
    %
    % Solve a 2D-Poisson equation at [0, 1] x [0, 1]
    % 
    %                      - u_xx - u_yy = f(x, y) 
    % subject to
    %                     u = 0 along the boundary
    % where
    % 
    %               f(x, y) = 2*pi*pi*sin(pi*x)*sin(pi*y)
    % 
    % Analytical solution:
    % 
    %               u(x, y) = g(x, y) = sin(pi*x)*sin(pi*y)
    %
    %   Discretization results in a linear system Au = b

    % Source term and boundary values
    f = @(x, y) 2 * pi * pi * sin( pi*x ) * sin( pi*y );
    g = @(x, y) sin( pi*x ) * sin( pi*y );
    NODES1D = 5;
    INNODES1D = NODES1D - 2;
    aStart = 0;
    aEnd = 1;
    h = ( aEnd - aStart ) / (NODES1D - 1); %;

    % Left-hand side (LHS) matrix 'A'
    % LHS is composed by the finite-difference coefficients for u''(x)
    M = 4*eye(INNODES1D) - diag(ones(INNODES1D-1,1), 1) ...
        -diag(ones(INNODES1D-1, 1), -1);
    N = -1*eye(INNODES1D);
    Z = zeros(INNODES1D);
    A = [M N Z; N M N; Z N M];
    % For INNODES1D > 3, previous line requires an external 'blktridiag' script
    % A = blktridiag(M, N, N, INNODES1D);

    % Right-hand side (RHS) 'b'
    % RHS is composed by the contributions from 
    % the boundary conditions g(x) and the source term f(x).
    b = zeros(INNODES1D*INNODES1D, 1);
    hSquare = h*h;
    for i = 1:INNODES1D
        for j = 1:INNODES1D
            b( (i-1)*INNODES1D + j, 1 ) = hSquare * f(i*h, j*h);
        end
    end

    % Exact solution 'uExact'
    uExact = ...
    [0.526514643772757;
   0.744604150011472;
   0.526514643772757;
   0.744604150011472;
   1.053029287545515;
   0.744604150011472;
   0.526514643772757;
   0.744604150011472;
   0.526514643772757];

    % CALL ALGORITHMS: PD_GMRES, ADAPTIVE_GMRES, SWITCH_GMRES, etc.
    
    % PD_GMRES
    tol = 1e-9;
    maxit = 100;
    mInitial = 3;
    uPD_GMRES = pd_gmres(A, b, mInitial, [], [], tol, maxit);
    
    % Assert whether the pd_gmres solution match the exact knonw solution
    assertElementsAlmostEqual(uPD_GMRES, uExact)
end
