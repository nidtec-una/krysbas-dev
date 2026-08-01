function [x, flag, relresvec, kdvec, time] = ...
    slgmres_e(A, b, m, l, d, epsilonThreshold, tol, maxit, xInitial, eigstol)
    % SLGMRES-E algorithm
    %
    %   Description:
    %   ------------
    %
    %   SLGMRES-E(m, l, d) is a restarted GMRES variant that switches, cycle
    %   by cycle, between two augmentation strategies based on an observed
    %   convergence-slowdown signal:
    %
    %   (a) By default, it behaves like LGMRES(m, l): the restart subspace is
    %       augmented with up to l error approximation vectors from prior
    %       cycles.
    %
    %   (b) When slowdown is detected -- the relative residual ratio between
    %       two consecutive cycles exceeds 1 - epsilonThreshold -- it
    %       switches to GMRES-E(m, d) for one or more cycles: the subspace
    %       is instead augmented with d harmonic Ritz vectors computed from
    %       the cycle that triggered the detection. It switches back to
    %       LGMRES-style augmentation as soon as the ratio improves again.
    %
    %   The slowdown signal, epsilon in ||r_m^(j)|| / ||r_m^(j-1)|| = 1 -
    %   epsilon, is read directly off the Givens-rotated residual already
    %   computed for the least-squares solve every cycle -- no extra cost.
    %   Under this switching rule, the residual norm is provably
    %   non-increasing cycle to cycle (Theorem 2 of [1]).
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relresvec, kdvec, time] = ...
    %       slgmres_e(A, b, m, l, d, epsilonThreshold, tol, maxit, xInitial)
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
    %   m:          int, optional
    %               Restart parameter, fixed across all cycles.  If m == n,
    %               the built-in unrestarted gmres is used.  Default:
    %               min(n, 10).
    %
    %   l:          int, optional
    %               Number of error approximation vectors to append during
    %               LGMRES-style cycles.  Must satisfy l > 0.  Default: 3.
    %
    %   d:          int, optional
    %               Number of harmonic Ritz vectors to append during
    %               GMRES-E-style cycles.  Must satisfy d > 0.  Default:
    %               min(m, 3).
    %
    %   epsilonThreshold:
    %               float, optional
    %               Slowdown threshold.  A cycle is classified as stagnating
    %               when ||r^(j)|| / ||r^(j-1)|| >= 1 - epsilonThreshold.
    %               Default: 0.01, as used in the numerical experiments of
    %               [1].
    %
    %   tol:        float, optional
    %               Relative residual tolerance.  Default: 1e-6.
    %
    %   maxit:      int, optional
    %               Maximum number of restart cycles.  Default: min(n, 10).
    %
    %   xInitial:   n-by-1 vector, optional
    %               Initial guess.  Default: zeros(n, 1).
    %
    %   eigstol:    float, optional
    %               Tolerance for the built-in eigs solver used inside
    %               harmonic_ritz_vectors.  Default: 1e-6.
    %
    %
    %   Output Parameters:
    %   ------------------
    %
    %   x:          n-by-1 vector
    %               Approximate solution to Ax = b.
    %
    %   flag:       integer (0 or 1)
    %               1 if the relative residual dropped below tol within
    %               maxit cycles; 0 otherwise.
    %
    %   relresvec:  (cycles+1)-by-1 vector
    %               Relative residual norm after each cycle, starting from 1.
    %
    %   kdvec:      (cycles+1)-by-1 vector
    %               Krylov subspace dimension used at each cycle (m during
    %               cycle 1, m+l or m+d thereafter depending on which
    %               augmentation was active).
    %
    %   time:       float
    %               Wall-clock time in seconds.
    %
    %
    %   References:
    %   -----------
    %
    %   [1] Cabral, J. C., Schaerer, C. E., & Bhaya, A. (2020). Improving
    %   GMRES(m) using an adaptive switching controller. Numerical Linear
    %   Algebra with Applications, 27(5), e2305.
    %
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

    % =========================================================================
    % ----> Sanity check on the number of input parameters
    % =========================================================================

    if nargin < 2
        error("Too few input parameters. Expected at least A and b.");
    elseif nargin > 10
        error("Too many input parameters.");
    end

    % =========================================================================
    % ----> Sanity checks on matrix A
    % =========================================================================

    if isempty(A)
        error("Matrix A cannot be empty.");
    end

    [rowsA, colsA] = size(A);
    if rowsA ~= colsA
        error("Matrix A must be square.");
    end

    n = rowsA;
    clear rowsA colsA;

    % =========================================================================
    % ----> Sanity checks on vector b
    % =========================================================================

    if isempty(b)
        error("Vector b cannot be empty.");
    end

    [rowsb, colsb] = size(b);
    if colsb ~= 1
        error("Vector b must be a column vector.");
    end

    if rowsb ~= n
        error("Dimension mismatch between matrix A and vector b.");
    end

    clear rowsb colsb;

    % =========================================================================
    % ----> Default value and sanity checks for m
    % =========================================================================

    if (nargin < 3) || isempty(m)
        m = min(n, 10);
    end

    if m == n
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b);
        time = toc();
        x = gmres_x;
        flag = (gmres_flag == 0);
        relresvec = resvec ./ resvec(1);
        kdvec = m .* ones(length(relresvec), 1);
        return
    end

    if m > n
        error("m must satisfy: 1 <= m <= n.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for l and d
    % =========================================================================

    if (nargin < 4) || isempty(l)
        l = 3;
    end

    if (nargin < 5) || isempty(d)
        d = min(m, 3);
    end

    if l <= 0
        error("l must satisfy: l > 0.");
    end

    if d <= 0
        error("d must satisfy: d > 0.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for epsilonThreshold
    % =========================================================================

    if (nargin < 6) || isempty(epsilonThreshold)
        epsilonThreshold = 0.01;
    end

    if epsilonThreshold <= 0 || epsilonThreshold >= 1
        error("epsilonThreshold must satisfy: 0 < epsilonThreshold < 1.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for tol
    % =========================================================================

    if (nargin < 7) || isempty(tol)
        tol = 1e-6;
    end

    if tol < eps
        warning("Tolerance is too small and it will be changed to eps.");
        tol = eps;
    elseif tol >= 1
        warning("Tolerance is too large and it will be changed to 1 - eps.");
        tol = 1 - eps;
    end

    % =========================================================================
    % ----> Default value for maxit
    % =========================================================================

    if (nargin < 8) || isempty(maxit)
        maxit = min(n, 10);
    end

    % =========================================================================
    % ----> Default value and sanity checks for xInitial
    % =========================================================================

    if (nargin < 9) || isempty(xInitial)
        xInitial = zeros(n, 1);
    end

    [rowsxInitial, colsxInitial] = size(xInitial);
    if colsxInitial ~= 1
        error("Initial guess xInitial is not a column vector.");
    end

    if rowsxInitial ~= n
        msg = "Dimension mismatch between matrix A and initial guess xInitial.";
        error(msg);
    end

    clear rowsxInitial colsxInitial;

    % =========================================================================
    % ----> Default value for eigstol
    % =========================================================================

    if (nargin < 10) || isempty(eigstol)
        eigstol = 1e-6;
    end

    % =========================================================================
    % ----> SLGMRES-E algorithm starts here
    % =========================================================================

    normY = 1 - epsilonThreshold; % stagnation threshold on the residual ratio

    flag = 0;
    x = xInitial;
    r0 = b - A * x;
    res1 = norm(r0); % initial residual norm; fixed denominator for relresvec

    relresvec = zeros(maxit + 1, 1);
    relresvec(1) = 1.0;
    kdvec = zeros(maxit + 1, 1);
    nCycles = 0;

    % Sliding window of LGMRES-style error approximation vectors, newest
    % vector last.  nZ tracks how many columns are actually populated
    % (independent of how many cycles have elapsed, since GMRES-E-style
    % cycles do not touch this window).
    zMat = zeros(n, l);
    nZ = 0;

    tic(); % start wall-clock timer

    % -------------------------------------------------------------------
    % Cycle 1: plain GMRES(m).  Neither augmentation strategy has any
    % history to draw on yet.
    % -------------------------------------------------------------------
    v1 = r0 / res1;
    [H, V, s] = modified_gram_schmidt_arnoldi(A, v1, m);
    [HUpTri, g] = plane_rotations(H, res1);

    Rs = HUpTri(1:s, 1:s);
    gs = g(1:s);
    minimizer = Rs \ gs;
    zCycle = V * minimizer;
    x = x + zCycle;
    nCycles = nCycles + 1;
    relresvec(nCycles + 1) = abs(g(s + 1)) / res1;
    kdvec(nCycles + 1) = s;

    if relresvec(nCycles + 1) < tol
        flag = 1;
        relresvec = relresvec(1:nCycles + 1);
        kdvec = kdvec(1:nCycles + 1);
        time = toc();
        return
    end

    % The correction from this cycle is always kept as the first error
    % approximation vector, regardless of what the next cycle turns out
    % to be.
    zMat(:, 1) = zCycle;
    nZ = 1;

    stagnating = (relresvec(nCycles + 1) / relresvec(nCycles) >= normY);
    if stagnating
        % Cycle 1's basis is a plain (non-augmented) Arnoldi basis, so the
        % cheap H'-based formula for Fold is exact here (equivalent to,
        % but cheaper than, W'*A'*W).  Later cycles cannot use this
        % shortcut because V's augmented columns get overwritten with raw
        % (non-Arnoldi-consistent) vectors below -- see the main loop.
        Fold = H(1:s, 1:s)';
        G = Rs' * Rs;
        dy = harmonic_ritz_vectors(Fold, G, d, V, eigstol);
    end

    % =========================================================================
    % ----> Main loop: cycles 2, 3, ...
    % =========================================================================

    while flag == 0 && nCycles < maxit

        r = b - A * x;
        beta = norm(r);
        v1 = r / beta;

        if ~stagnating
            % --- LGMRES-style cycle ---
            % Note the fliplr placement here is the opposite of the
            % GMRES-E-style branch below: matching lgmres.m's own
            % convention (not gmres_e.m's), the newest error vector goes
            % into the first augmentation slot V(:,m+1).
            lUse = min(nZ, l);
            [H, V, s] = ...
                augmented_gram_schmidt_arnoldi(A, v1, m, zMat(:, 1:lUse));
            [HUpTri, g] = plane_rotations(H, beta);
            Rs = HUpTri(1:s, 1:s);
            gs = g(1:s);
            minimizer = Rs \ gs;
            V(:, m + 1:s) = fliplr(zMat(:, 1:lUse));
        else
            % --- GMRES-E-style cycle ---
            [H, V, s] = ...
                augmented_gram_schmidt_arnoldi(A, v1, m, fliplr(dy(:, 1:d)));
            [HUpTri, g] = plane_rotations(H, beta);
            Rs = HUpTri(1:s, 1:s);
            gs = g(1:s);
            minimizer = Rs \ gs;
            V(:, m + 1:s) = dy(:, 1:d);
        end

        aux = V * minimizer;
        x = x + aux;
        nCycles = nCycles + 1;
        relresvec(nCycles + 1) = abs(g(s + 1)) / res1;
        kdvec(nCycles + 1) = s;

        if relresvec(nCycles + 1) < tol
            flag = 1;
            break
        end

        % ------------------------------------------------------------------
        % Decide the augmentation strategy for the NEXT cycle from the
        % ratio just observed.  This mirrors [1] eq. (33)-(35): the ratio
        % is read off the residual we already computed above, at no extra
        % cost.  The decision is independent of which strategy produced
        % THIS cycle's result -- a recovering GMRES-E cycle's correction
        % is just as valid an error-approximation vector as an LGMRES
        % cycle's, and a newly-stagnating LGMRES cycle's basis is just as
        % valid a source of harmonic Ritz vectors as a GMRES-E cycle's.
        % ------------------------------------------------------------------
        ratio = relresvec(nCycles + 1) / relresvec(nCycles);

        if ratio >= normY
            stagnating = true;
            W = V(:, 1:s);
            Fold = W' * A' * W;
            G = Rs' * Rs;
            dy = harmonic_ritz_vectors(Fold, G, d, V, eigstol);
        else
            stagnating = false;
            if nZ < l
                nZ = nZ + 1;
                zMat(:, nZ) = aux;
            else
                zMat(:, 1:l - 1) = zMat(:, 2:l);
                zMat(:, l) = aux;
            end
        end

    end

    relresvec = relresvec(1:nCycles + 1);
    kdvec = kdvec(1:nCycles + 1);
    time = toc();

end
