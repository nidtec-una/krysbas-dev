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

function test_poisson_one_dimension_dirichlet_bc()
    % Test solvers for a 1D-FDM linear system with Dirichlet boundaries.
    %
    % Description:
    % ============
    %   Solve a 1D-Poisson equation in the domain [aStart, aEnd] = [0, 1]
    %
    %                         - u''(x) = f(x)
    %
    %   subject to u(0) = g(0) = 0, u(1) = g(1) = 1; where f(x) = 6*x and
    %   g(x) = x^3.
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
    % b(1, 1) = g( astart );
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
