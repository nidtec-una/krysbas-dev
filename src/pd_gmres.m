function [x, flag, relres, iter, resvec, restarted, time] = ...
    pd_gmres(A, b, mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD, ...
             varargin)
    % PD-GMRES Proportional-Derivative GMRES(m)
    %
    %   pd_gmres is a modified implementation of the restarted Generalized
    %   Minimal Residual Error or GMRES(m) [1], performed by using a
    %   proportional-derivative control-inspired law to update adaptively the
    %   restarting parameter m before each restart. This implementation follows
    %   closely the one presented in [2].
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relres, iter, resvec, restarted, time] = pd_gmres(A, b, ...
    %       mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD)
    %
    %
    %   Input Parameters:
    %   -----------------
    %
    %   A:        n-by-n matrix
    %             Left-hand side of the linear system Ax = b.
    %
    %   b:        n-by-1 vector
    %             Right-hand side of the linear system Ax = b.
    %
    %   mInitial: int, optional
    %             Initial restart parameter (similar to 'restart' in MATLAB). If
    %             empty, not given, or equal to n, then mInitial is not used and
    %             the full unrestarted gmres algorithm with maxit = min(n, 10)
    %             is employed. Note that we require 1 <= mInitial <= n.
    %
    %   mMinMax:  2-by-1 vector, optional
    %             Minimum and maximum values of the restart paramter m. Default
    %             is [1; n]. We require 1 <= mMinMax(1) < mMinMax(2) <= n. Note
    %             that for large matrices (e.g., > 600), the default value
    %             mMinMax(2) = n might require large computational resources.
    %
    %   mStep:    int, optional
    %             Step size for increasing the restart parameter m when
    %             m < mMinMax(1). Default is 1 if n <= 10 and 3 otherwise.
    %
    %   tol:      float, optional
    %             Tolerance error threshold for the relative residual norm.
    %             Default is 1e-6.
    %
    %   maxit:    int, optional
    %             Maximum number of outer iterations.
    %
    %   xInitial: n-by-1 vector, optional
    %             Vector of initial guess. Default is zeros(n, 1).
    %
    %   alphaPD:  2-by-1 vector, optional
    %             Proportional and derivative coefficients from the
    %             proportional-derivative controller. Default is [-3, 5],
    %             following Ref. [2].
    %
    %   Output parameters:
    %   ------------------
    %
    %   x:        n-by-1 vector
    %             Approximate solution of the linear system.
    %
    %   flag:     boolean
    %             1 if the algorithm converged, 0 otherwise.
    %
    %
    %   log_res:  (1 up to to max_iter)-by-1 vector
    %             relative residual norms
    %
    %   References:
    %   -----------
    %
    %   [1] Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal
    %   residual algorithm for solving nonsymmetric linear systems. SIAM
    %   Journal on scientific and statistical computing, 7(3), 856-869.
    %
    %   [2] Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A
    %   proportional-derivative control strategy for restarting the GMRES(m)
    %   algorithm. Journal of Computational and Applied Mathematics,
    %   337, 209-224.
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

    % ----> Sanity check on the number of input parameters
    if nargin < 2
        error("Too few input parameters. Expected at least A and b.");
    elseif nargin > 9
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

    % ----> Default values and sanity checks for parameter mInitial
    % When the restart parameter is not specified explicitly as an input
    % argument, pd_gmres() behaves identically to the unrestarted gmres()
    % with a fixed number of iterations given by min(n, 10).
    %
    % There are three possibilities when this might happen:
    %   (1) If only two arguments are given i.e., A and b.
    %   (2) If mInitial equals the dimension of A i.e., mInitial = n.
    %   (3) If an empty matrix is given as mInitial, i.e., mInitial = [].
    if (nargin < 3) || isempty(mInitial) || (mInitial == n)
        restarted = false;  % use unrestarted pd_gmres()
    else
        restarted = true;  % use restarted pd_gmres()
    end

    % If the restarted version of pd_gmres will be used, then the value of
    % mInitial must be bounded between 1 and n.
    if restarted && (mInitial < 1 || mInitial > n)
        error("mInitial must satisfy: 1 <= mInitial <= n.");
    end

    % ----> Default values and sanity checks for mMinMax
    if (nargin < 4) || isempty(mMinMax)
        mMinMax = [1; n];
    else
        if ~restarted
            warning("mMinMax was given but will not be used.");
        end
    end

    if restarted && ...
            ((mMinMax(1) < 1) || ...
             (mMinMax(2) > n) || ...
             (mMinMax(2) <= mMinMax(1)))
        error("mMinMax must satisfy: 1 <= mMinMax(1) < mMinMax(2) <= n.");
    end

    if restarted && ((mMinMax(1) > mInitial) || (mMinMax(2) < mInitial))
        error("mMinMax must satisfy: mMinMax(1) <= mInitial <= mMinMax(2).");
    end

    mMin = mMinMax(1);
    mMax = mMinMax(2);

    % ----> Default values and sanity checks for mStep
    if (nargin < 5) || isempty(mStep)
        mStep = 1;
    else
        if ~restarted
            warning("mStep was given but will not be used.");
        end
    end

    if (mStep < 1) || (mStep > n - 1)
        error("mStep must satisfy: 0 < mStep < n.");
    end

    % ----> Default value and sanity checks for tol
    if (nargin < 6) || isempty(tol)
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
    if (nargin < 7) || isempty(maxit)
        if restarted
            maxit = min(ceil(n / mInitial), 10);
        else
            maxit = min(n, 10);
        end
    end

    % ----> Default value and sanity checks for initial guess xInitial
    if (nargin < 8) || isempty(xInitial)
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

    % ----> Default value for alphaPD
    % The default values for these parameteres were taken from page 217 of [2]
    if (nargin < 9) || isempty(alphaPD)
        alphaPD = [-3; 5];
    end

    alphaP = alphaPD(1);
    alphaD = alphaPD(2);

    % ---> Algorithm starts here

    % Unrestarted version of the PD-GMRES
    % This block calls the bult-in gmres() function but outputs slightly
    % different variables. In particular, the resvec vector only shows the
    % initial and the last residual (contrary to gmres() that stores all the
    % residuals).
    if ~restarted
        tic();
        [x, flag, relres, iter, resvec] = ...
            gmres(A, b, [], tol, maxit, [], [], xInitial);
        resvec = [resvec(1); resvec(end)];
        time = toc();
        return
    end

    % Restarted version of PD-GMRES

    % Algorithm setup
    flag = 0;
    restart = 1;
    r0 = b - A * xInitial;
    res(1, :) = norm(r0);
    resvec(1, :) = (norm(r0) / res(1, 1));
    iter(1, :) = restart;
    mIteration(1, 1) = mInitial;

    tic();  % start measuring CPU time

    while flag == 0

        if iter(size(iter, 1), :) ~= 1
            [miter] = pd_rule(m, mInitial, mMin, res, ...
                              iter(size(iter, 1), :), ...
                              mStep, mMax, alphaP, alphaD);
            m = miter(1, 1);
            mInitial = miter(1, 2);
        else
            m = mInitial;
        end

        mIteration(iter(size(iter, 1), :) + 1, 1) = m;
        r = b - A * xInitial;
        beta = norm(r);
        v1 = r / beta;

        % Apply modified Gram-Schmidt Arnoldi
        [H, V] = modified_gram_schmidt_arnoldi(A, v1, m);

        g = zeros(m + 1, 1);
        g(1, 1) = beta;

        % Plane rotations (QR decompostion)
        for j = 1:m
            P = eye(m + 1);
            sin = H(j + 1, j) / (sqrt(H(j + 1, j)^2 + H(j, j)^2));
            cos = H(j, j) / (sqrt(H(j + 1, j)^2 + H(j, j)^2));
            P(j, j) = cos;
            P(j + 1, j + 1) = cos;
            P(j, j + 1) = sin;
            P(j + 1, j) = -sin;
            H = P * H;
            g = P * g;
        end

        R = zeros(m, m);
        G = zeros(m, 1);

        for k = 1:m
            G(k) = g(k);
            for i = 1:m
                R(k, i) = H(k, i);
            end
        end

        minimizer = R \ G;
        xm = xInitial + V(1:n, 1:m) * minimizer;
        res(restart + 1, :) = abs(g(m + 1, 1));
        iter(restart + 1, :) = restart + 1;
        resvec(size(resvec, 1) + 1, :) = abs(g(m + 1, 1) / res(1, 1));
        % Use last component of g as residual
        if abs(g(m + 1, 1)) / res(1, 1) < tol || size(resvec, 1) == maxit
            flag = 1;  % solution has converged
            x = xm;  % solution vector
        else
            xInitial = xm;  % update and restart
            restart = restart + 1;
        end
        % Compute the relative residual
        relres = norm(b - A * xm) / norm(b);
    end

    time = toc();  % record CPU time

end
