function [x, flag, relresvec, kdvec, time] = ...
    gmres_e(A, b, m, d, tol, maxit, xInitial, eigstol, varargin)
    % GMRES-E algorithm
    %
    %   GMRES-E is a modified implementation of the restarted
    %   Generalized Minimal Residual Error or GMRES(m) [1], performed by
    %   appending 'd' eigenvectors corresponding to a few of the smallest
    %   eigenvalues in magnitude for each outer iteration. In practice, the
    %   approximate eigenvectors are the harmonic Ritz vectors associated to
    %   the harmonic Ritz values per outer iteration.
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relresvec, time] = ...
    %       gmres_e(A, b, m, k, tol, maxit, xInitial, eigstol)
    %
    %
    %   Input Parameters:
    %   -----------------
    %
    %   A:          n-by-n matrix
    %               Left-hand side of the linear system Ax = b.
    %
    %   b:          n-by-1 vector
    %               Right-hand side of the linear system Ax = b.
    %
    %   m:          int
    %               Restart parameter (similar to 'restart' in MATLAB).
    %               If m == n, built-in unrestarted gmres will be used
    %
    %   d:          int
    %               Number of eigenvectors corresponding to a few of the
    %               smallest eigenvalues in magnitude for each outer
    %               iteration. Default is min(m, 3), but values between
    %               1 and 5 are typical. According to [1], "even just a few
    %               eigenvectors can make a big difference if the matrix has
    %               both small and large eigenvalues". If m < n AND d == 0,
    %               the built-in gmres(m) will be used. If d > m, an error
    %               is raised.
    %
    %   tol:        float, optional
    %               Tolerance error threshold for the relative residual norm.
    %               Default is 1e-6.
    %
    %   maxit:      int, optional
    %               Maximum number of outer iterations.
    %               Default is min(m, 10).
    %
    %   xInitial:   n-by-1 vector, optional
    %               Vector of initial guess. Default is zeros(n, 1).
    %
    %   eigstol:    float, optional
    %               Tolerance for computing eigenvectors using the built-in
    %               MATLAB function `eigs`. Default is 1e-6.
    %
    %
    %   Output parameters:
    %   ------------------
    %
    %   x:          n-by-1 vector
    %               Approximate solution to the linear system.
    %
    %   flag:       boolean
    %               1 if the algorithm has converged, 0 otherwise.
    %
    %   relresvec:  (1 up to maxit)-by-1 vector
    %               Vector of relative residual norms of every outer
    %               iteration (cycles). The last relative residual norm is
    %               simply given by relresvec(end).
    %
    %   kdvec:      (1 up to maxit)-by-1 vector
    %               For LGMRES, kdvec is a constant vector whose elements
    %               correspond to the size of the Krylov subspace, i.e.,
    %               m + d. Note that in some cases, there could be a
    %               "happy breakdown" where the dimension of the search
    %               space < m + d.
    %
    %   time:       float
    %               Computational time in seconds.
    %
    %   References:
    %   -----------
    %
    %   [1] Morgan, R. B. (1995). A restarted GMRES method augmented with
    %   eigenvectors. SIAM Journal on Matrix Analysis and Applications,
    %   16(4), 1154-1171.
    %
    %   Copyright:
    %   ----------
    %
    %   This file is part of the KrySBAS MATLAB Toolbox.
    %
    %   Copyright 2023 CC&MA - NIDTec - FP - UNA
    %
    %   KrySBAS is free software: you can redistribute it and/or modify it
    %   under the terms of the GNU General Public License as published by the
    %   Free Software Foundation, either version 3 of the License, or (at
    %   your option) any later version.
    %
    %   KrySBAS is distributed in the hope that it will be useful, but WITHOUT
    %   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    %   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    %   for more details.
    %
    %   You should have received a copy of the GNU General Public License along
    %   with this file.  If not, see <http://www.gnu.org/licenses/>.
    %

    % ----> Sanity check on the number of input parameters
    if nargin < 2
        error("Too few input parameters. Expected at least A and b.");
    elseif nargin > 8
        error("Too many input parameters.");
    end

    % ----> Sanity checks on matrix A
    % Check whether A is non-empty
    if isempty(A)
        error("Matrix A cannot be empty.");
    end

    % Check whether A is square
    [rowsA, colsA] = size(A);
    if rowsA ~= colsA
        error("Matrix A must be square.");
    end

    n = rowsA;
    clear rowsA colsA;

    % ----> Sanity checks on vector b
    % Check whether b is non-empty
    if isempty(b)
        error("Vector b cannot be empty.");
    end

    % Check whether b is a column vector
    [rowsb, colsb] = size(b);
    if colsb ~= 1
        error("Vector b must be a column vector.");
    end

    % Check whether b has the same number of rows as b
    if rowsb ~= n
        error("Dimension mismatch between matrix A and vector b.");
    end

    clear rowsb colsb;

    % Special sanity checks for GMRES-E here

    % ----> Default value and sanity checks for m
    if (nargin < 3) || isempty(m)
        m = min(n, 10);
    end

    % ----> If m > n, error message is printed
    if m > n
        error("m must satisfy: 1 <= m <= n.");
    end

    % ----> If m == n, built-in unrestarted gmres will be used
    if m == n
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b);
        time = toc();
        x = gmres_x;
        if gmres_flag == 0
            flag = 1;
        else
            flag = 0;
        end
        relresvec = resvec ./ resvec(1, 1);
        kdvec = m .* ones(length(relresvec), 1);
        return
    end

    % ----> If m < n AND d == 0, built-in gmres(m) will be used
    if (m < n) && (d == 0)
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b, m);
        time = toc();
        x = gmres_x;
        if gmres_flag == 0
            flag = 1;
        else
            flag = 0;
        end
        relresvec = resvec ./ resvec(1, 1);
        kdvec = m .* ones(length(relresvec), 1);
        return
    end

    % ----> Default value and sanity checks for d
    if (nargin < 4) || isempty(d)
        d = min(m, 3);
    end

    % ----> d cannot be larger than m
    if d > m
        error("d cannot be larger than m");
    end

    % Default value and sanity checks for tol
    if (nargin < 5) || isempty(tol)
        tol = 1e-6;
    end

    if tol < eps
        warning("Tolerance is too small and it will be changed to eps.");
        tol = eps;
    elseif tol >= 1
        warning("Tolerance is too large and it will be changed to 1-eps.");
        tol = 1 - eps;
    end

    % ----> Default value for maxit
    if (nargin < 6) || isempty(maxit)
        maxit = min(n, 10);
    end

    % ----> Default value and sanity checks for initial guess xInitial
    if (nargin < 7) || isempty(xInitial)
        xInitial = zeros(n, 1);
    end

    % Check whether xInitial is a column vector
    [rowsxInitial, colsxInitial] = size(xInitial);
    if colsxInitial ~= 1
        error("Initial guess xInitial is not a column vector.");
    end

    % Check whether x0 has the right dimension
    if rowsxInitial ~= n
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        error(msg);
    end

    clear rowsxInitial colsxInitial;

    % Default value for eigstol
    if (nargin < 8) || isempty(eigstol)
        eigstol = 1e-6;
    end

    % ---> GMRES-E algorithm starts here
    % First outer iteration is a simple restarted GMRES(m + k) execution
    % Afterwards, we can compute the 'd' eigenvectors if convergence is
    % not achieved.

    % Compute normalized residual vector
    flag = 0;
    restart = 1;
    r0 = b - A * xInitial;
    res(1, :) = norm(r0);
    relresvec(1, :) = (norm(r0) / res(1, 1));
    iter(1, :) = restart;
    beta = norm(r0);
    v1 = r0 / beta;

    tic(); % start measuring CPU time

    % Modified Gram-Schmidt Arnoldi iteration
    % This is the first run. Since we don't have
    % harmonic Ritz vectors yet, we use GMRES(s) = GMRES(m+d).
    s = m + d;
    [H, V, sUpdated] = modified_gram_schmidt_arnoldi(A, v1, s);

    % Plane rotations
    [HUpTri, g] = plane_rotations(H, beta);

    % Solve the least-squares problem
    s = sUpdated;
    Rs = HUpTri(1:s, 1:s);
    gs = g(1:s);
    minimizer = Rs \ gs;
    xm = xInitial + V * minimizer;

    % Update residual norm, iterations, and relative residual vector
    res(restart + 1, :) = abs(g(s + 1, 1));
    iter(restart + 1, :) = restart + 1;
    relresvec(restart + 1, :) = res(restart + 1, :) / res(1, 1);

    % Check convergence
    if relresvec(restart + 1, :) < tol
        % We reached convergence.
        flag = 1;
        x = xm;
        time = toc();
        kdvec = s .* ones(length(relresvec), 1);
        return
    else
        % We have not reached convergence.
        % For GMRES-E, an eigenvalue problem is solved.
        % Eigenvalue problem setup, from [1], p. 1161, step 5
        Fold = H(1:s, 1:s)';
        G = Rs' * Rs;
        dy = harmonic_ritz_vectors(Fold, G, d, V, eigstol);

        % Update and restart.
        restart = restart + 1;
    end

    % Empty matrix E for next outer iteration
    E = zeros(s, d);

    % ---> GMRES-E Algorithm for restart > 1

    while flag == 0 && restart <= maxit

        % Compute normalized residual vector
        r = b - A * xm;
        beta = norm(r);
        v1 = r / beta;

        % Augmented Gram-Schmidt Arnoldi iteration
        % Notice that for GMRES-E, we need E,
        % a n-by-k matrix, the matrix of eigenvectors
        % from the last outer iteration.
        [H, V, s] = ...
            augmented_gram_schmidt_arnoldi(A, v1, m, fliplr(dy(:, 1:d)));
        % Patch: to be solved
        % Why does the order matter?

        % Plane rotations
        [HUpTri, g] = plane_rotations(H, beta);

        % Solve the least-squares problem
        Rs = HUpTri(1:s, 1:s);
        gs = g(1:s);
        minimizer = Rs \ gs;

        % Replace last k vectors from matrix V with the approximate
        % eigenvectors, and compute the new approximate solution, as done
        % in step 4, p. 1161 of [1].
        V(:, m + 1:s) = dy(:, 1:d);
        x = xm + V * minimizer;

        % Update residual norm, iterations, and relative residual vector
        res(restart + 1, :) = abs(g(s + 1, 1));
        iter(restart + 1, :) = restart + 1;
        relresvec(restart + 1, :) = ...
            res(restart + 1, :) / res(1, 1);

        % Check convergence
        if relresvec(restart + 1, 1) < tol
            % We reached convergence.
            flag = 1;
            time = toc();
            kdvec = s .* ones(length(relresvec), 1);
            return

        elseif restart < maxit
            % We have not reached convergence.
            % Eigenvalue problem setup, from [1], p. 1161, step 5
            W = V(1:n, 1:s);
            Fold = W' * A' * W;
            G = Rs' * Rs;
            dy = harmonic_ritz_vectors(Fold, G, d, V, tol);
        end

        % Update and restart.
        xm = x;
        restart = restart + 1;
    end
    kdvec = s .* ones(length(relresvec), 1);
    time = toc();
end
