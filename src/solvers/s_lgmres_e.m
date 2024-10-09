function [x, flag, relresvec, kdvec, time] = ...
    s_lgmres_e(A, b, mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD, ...
             eps0, l, d, varargin)
    % Adaptive SLGMRES-E with Proportional-Derivative Controller
    %
    %   Description:
    %   ------------
    %
    %   Adaptive SLGMRES-E a hybrid method for improving GMRES(m) by
    %   restraining the stagnation and slowdown convergence. The idea behind
    %   the variation of the structure of the GMRES(m) relies on the fact that
    %   once either the slowdown of convergence or stagnation is detected,
    %   a rule of switching between two strategies chooses the one most
    %   appropriate for facing the identified convergence problem. For this
    %   setup. two monotonically decreasing strategies, obtained from the
    %   literature, are combined: the LGMRES and the GMRES-E, this results
    %   in the switching method denoted as SLGMRES-E. In some cases, a
    %   proportional derivative rule is introduced to modify the restart
    %   parameter 'm', which gives place to this Adaptive SLGMRES-E
    %   algorithm.
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relresvec, kdvec, time] = s_lgmres_e(A, b, ...
    %       mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD)
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
    %   mInitial:   int, optional
    %               Initial restart parameter (similar to 'restart' in MATLAB).
    %               If empty, not given, or equal to n, then mInitial is not
    %               used and the full unrestarted gmres algorithm with
    %               maxit = min(n, 10) is employed. Note that we require
    %               1 <= mInitial <= n.
    %
    %   mMinMax:    2-by-1 vector, optional
    %               Minimum and maximum values of the restart paramter m.
    %               Default is [1; n]. We require
    %               1 <= mMinMax(1) < mMinMax(2) <= n. Note that for large
    %               matrices (e.g., > 600), the default value mMinMax(2) = n
    %               might require large computational resources.
    %
    %   mStep:      int, optional
    %               Step size for increasing the restart parameter m when
    %               m < mMinMax(1). Default is 1 if n <= 10 and 3 otherwise.
    %
    %   tol:        float, optional
    %               Tolerance error threshold for the relative residual norm.
    %               Default is 1e-6.
    %
    %   maxit:      int, optional
    %               Maximum number of outer iterations.
    %
    %   xInitial:   n-by-1 vector, optional
    %               Vector of initial guess. Default is zeros(n, 1).
    %
    %   alphaPD:    2-by-1 vector, optional
    %               Proportional and derivative coefficients from the
    %               proportional-derivative controller. Default is [-3, 5],
    %               following Ref. [2].
    %
    %   eps0:       float
    %               Stagnation threshold parameter, as defined in [1].
    %               Slow convergence occurs when ||y|| < eps0,
    %               and approximate error vectors don't have enough
    %               impact.
    %
    %   Output parameters:
    %   ------------------
    %
    %   x:          n-by-1 vector
    %               Approximate solution of the linear system.
    %
    %   flag:       boolean
    %               1 if the algorithm has converged, 0 otherwise.
    %
    %   relressvec: (1 up to maxit)-by-1 vector
    %               Vector of relative residual norms of every outer iteration
    %               (cycles). The last relative residual norm is simply given
    %               by relresvec(end).
    %
    %   kdvec:      (1 up to maxit)-by-1 vector
    %               Vector of restart parameter values. In case the unrestarted
    %               algorithm is invoked, kdvec = NaN.
    %
    %   time:       scalar
    %               Computational time in seconds.
    %
    %   References:
    %   -----------
    %
    %   [1] Cabral, J. C., Schaerer, C. E., & Bhaya, A. (2020). Improving
    %   GMRES(m) using an adaptive switching controller. Numerical Linear
    %   Algebra with Applications, 27(5), e2305.
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
        restarted = false;
    else
        restarted = true;
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

    % ----> Default value for eps0
    if (nargin < 9) || isempty(eps0)
        eps0 = 0.01;
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

    % ---> S-LGMRES-E Algorithm starts here
    % First outer iteration consists of a single GMRES(m) iteration, plus a
    % switching strategy:
    %
    % * It it presents an appropriate good convergence, computing
    %   approximation error vectors is enough, as in the LGMRES mehtod
    % * If poor convergence is rather detected, then the Harmonic Ritz 
    %   vectors are computed as in the GMRES-E method, plus an update of the
    %   restart parameter 'm' via a proportional-derivative (PD) rule
    
    flag = 0;
    flag_stagnation = 0;
    restart = 1;
    r0 = b - A * xInitial;
    res(1, :) = norm(r0);
    relresvec(1, :) = (norm(r0) / res(1, 1));
    iter(1, :) = restart;
    beta = norm(r0);
    v1 = r0 / beta;
    kdvec(1, 1) = mInitial;

    tic();  % start measuring CPU time

    % Modified Gram-Schmidt Arnoldi iteration
    % This is the first run. Since we don't have
    % harmonic Ritz vectors yet, we use GMRES(m).
    s = mInitial;
    [H, V, sUpdated] = modified_gram_schmidt_arnoldi(A, v1, s);
    
    % Plane rotations
    [HUpTri, g] = plane_rotations(H, beta);

    % Solve the least-squares problem
    s = sUpdated;
    Rs = HUpTri(1:s, 1:s);
    gs = g(1:s);
    minimizer = Rs \ gs;
    zCurrentCycle = V * minimizer;
    xm = xInitial + zCurrentCycle;

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
        % Now, we check for stagnation of the 'minimizer' vector
        if norm(minimizer) < eps0
            flag_stagnation = 1;
            Fold = H(1:s, 1:s)';
            G = Rs' * Rs;
            dy = harmonic_ritz_vectors(Fold, G, d, V, eigstol);

            % compute pdrule
        else
            % save zCurrentcycle
        end

        Fold = H(1:s, 1:s)';
        G = Rs' * Rs;
        dy = harmonic_ritz_vectors(Fold, G, d, V, eigstol);

        % Control block for 'm'
        if iter(size(iter, 1), :) ~= 1
            [miter] = pd_rule(m, n, mInitial, mMin, mMax, ...
                              mStep, res, iter(size(iter, 1), :), ...
                              alphaP, alphaD);
            m = miter(1, 1);
            mInitial = miter(1, 2);
        else
            m = mInitial;
        end
        kdvec(iter(size(iter, 1), :) + 1, 1) = m;

        % Update and restart.
        restart = restart + 1;
    end

    % Empty matrix E for next outer iteration
    E = zeros(s, d);
    
    % ---> S-LGMRES-E Algorithm for restart > 1
    while flag == 0 && restart <= maxit

        % Compute normalized residual vector
        r = b - A * xm;
        beta = norm(r);
        v1 = r / beta;

        % Now we check for the stagnation (or not) of the previous cycle
        % * If flag_stagnation == 0, then we generate the Arnoldi basis using
        %   the error approximation vectors,
        % * If flag_stagnation == 1, then we generate the Arnoldi basis using
        %   the harmonic Ritz vectors.
        if flag_stagnation == 0
            %s = mj + l;
            [H, V, s] = ...
                augmented_gram_schmidt_arnoldi(A, v1, m, zMat(:, 1:min(nevec, l)));
        else
            %s = mj + d;
            [H, V, s] = ...
                augmented_gram_schmidt_arnoldi(A, v1, m, fliplr(dy(:, 1:d)));
        end

        % To be continued...

    end
    
    time = toc();

end
