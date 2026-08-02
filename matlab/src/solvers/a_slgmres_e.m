function [x, flag, relresvec, kdvec, time] = ...
    a_slgmres_e(A, b, mInitial, mMinMax, mStep, l, d, epsilonThreshold, ...
                alphaPD, tol, maxit, xInitial, eigstol, varargin)
    % A-SLGMRES-E algorithm
    %
    %   Description:
    %   ------------
    %
    %   A-SLGMRES-E(mj, l, d) extends slgmres_e (SLGMRES-E) with an adaptive
    %   restart parameter: whenever a cycle triggers the same
    %   convergence-slowdown signal that switches slgmres_e from LGMRES-style
    %   to GMRES-E-style augmentation, this variant ALSO grows the restart
    %   parameter m via the proportional-derivative law of pd_rule (the same
    %   law used by pd_gmres), for as long as slowdown persists.  m only
    %   changes on stagnating cycles; it stays fixed during LGMRES-style
    %   cycles.  This is Algorithm 1 of [1].
    %
    %   See slgmres_e for the switching rule itself (eq. 33-35 of [1]) and
    %   why it can be read off the Givens-rotated residual at no extra cost.
    %
    %   Deviation from the reference implementation: pd_rule's "warm-up"
    %   gating (it only applies the derivative term once at least 3 cycles
    %   of residual history exist, and only the proportional term with 2)
    %   is driven here by the complete, correctly-indexed cycle history.
    %   Cabral's own reference script (Adaptive_PD_lgmres_e.m) has a
    %   residual-history update for its first cycle commented out, which
    %   leaves its own cycle counter one cycle behind for the rest of the
    %   run -- traced and confirmed by reproducing that exact lag in a
    %   scratch copy of this port, which closed most (105 vs. 102 cycles,
    %   from 118) but not all of the gap to the reference on the sherman5
    %   case used in [1]'s own numerical experiments. Both this port and
    %   the reference converge correctly; they simply grow m on a slightly
    %   different schedule.
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relresvec, kdvec, time] = ...
    %       a_slgmres_e(A, b, mInitial, mMinMax, mStep, l, d, ...
    %                   epsilonThreshold, alphaPD)
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
    %               Initial restart parameter.  If mInitial == n, the
    %               built-in unrestarted gmres is used.  Default: min(n, 10).
    %
    %   mMinMax:    2-by-1 vector, optional
    %               Minimum and maximum values the restart parameter m may
    %               take.  Default: [1; n].  We require
    %               1 <= mMinMax(1) < mMinMax(2) <= n.
    %
    %   mStep:      int, optional
    %               Step size for increasing the internal mInitial tracked
    %               by pd_rule when m falls below mMinMax(1).  Default: 1.
    %
    %   l:          int, optional
    %               Number of error approximation vectors to append during
    %               LGMRES-style cycles.  Must satisfy l > 0.  Default: 3.
    %
    %   d:          int, optional
    %               Number of harmonic Ritz vectors to append during
    %               GMRES-E-style cycles.  Must satisfy d > 0.  Default:
    %               min(mInitial, 3).
    %
    %   epsilonThreshold:
    %               float, optional
    %               Slowdown threshold.  A cycle is classified as stagnating
    %               when ||r^(j)|| / ||r^(j-1)|| >= 1 - epsilonThreshold.
    %               Default: 0.01, as used in the numerical experiments of
    %               [1].
    %
    %   alphaPD:    2-by-1 vector, optional
    %               Proportional and derivative coefficients for the
    %               restart-parameter growth law (pd_rule).  Default:
    %               [2; 0.8], the values reported in [1] -- note this
    %               differs from pd_gmres's own default of [-3; 5], which
    %               comes from a different paper ([2]) tuned for a different
    %               purpose (shrinking as well as growing m every cycle,
    %               rather than only growing it on detected stagnation).
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
    %               Krylov subspace dimension used at each cycle.
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
    %   [2] Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A
    %   proportional-derivative control strategy for restarting the GMRES(m)
    %   algorithm. Journal of Computational and Applied Mathematics,
    %   337, 209-224.
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
    elseif nargin > 13
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
    % ----> Default value and sanity checks for mInitial
    % =========================================================================

    if (nargin < 3) || isempty(mInitial)
        mInitial = min(n, 10);
    end

    if mInitial == n
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b);
        time = toc();
        x = gmres_x;
        flag = (gmres_flag == 0);
        relresvec = resvec ./ resvec(1);
        kdvec = mInitial .* ones(length(relresvec), 1);
        return
    end

    if mInitial < 1 || mInitial > n
        error("mInitial must satisfy: 1 <= mInitial <= n.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for mMinMax
    % =========================================================================

    if (nargin < 4) || isempty(mMinMax)
        mMinMax = [1; n];
    end

    if (mMinMax(1) < 1) || (mMinMax(2) > n) || (mMinMax(2) <= mMinMax(1))
        error("mMinMax must satisfy: 1 <= mMinMax(1) < mMinMax(2) <= n.");
    end

    if (mMinMax(1) > mInitial) || (mMinMax(2) < mInitial)
        error("mMinMax must satisfy: mMinMax(1) <= mInitial <= mMinMax(2).");
    end

    mMin = mMinMax(1);
    mMax = mMinMax(2);

    % =========================================================================
    % ----> Default value and sanity checks for mStep
    % =========================================================================

    if (nargin < 5) || isempty(mStep)
        mStep = 1;
    end

    if (mStep < 1) || (mStep > n - 1)
        error("mStep must satisfy: 0 < mStep < n.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for l and d
    % =========================================================================

    if (nargin < 6) || isempty(l)
        l = 3;
    end

    if (nargin < 7) || isempty(d)
        d = min(mInitial, 3);
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

    if (nargin < 8) || isempty(epsilonThreshold)
        epsilonThreshold = 0.01;
    end

    if epsilonThreshold <= 0 || epsilonThreshold >= 1
        error("epsilonThreshold must satisfy: 0 < epsilonThreshold < 1.");
    end

    % =========================================================================
    % ----> Default value for alphaPD
    % =========================================================================

    if (nargin < 9) || isempty(alphaPD)
        alphaPD = [2; 0.8];
    end

    alphaP = alphaPD(1);
    alphaD = alphaPD(2);

    % =========================================================================
    % ----> Default value and sanity checks for tol
    % =========================================================================

    if (nargin < 10) || isempty(tol)
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

    if (nargin < 11) || isempty(maxit)
        maxit = min(n, 10);
    end

    % =========================================================================
    % ----> Default value and sanity checks for xInitial
    % =========================================================================

    if (nargin < 12) || isempty(xInitial)
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

    if (nargin < 13) || isempty(eigstol)
        eigstol = 1e-6;
    end

    % =========================================================================
    % ----> A-SLGMRES-E algorithm starts here
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

    % Current restart parameter and pd_rule's own "initial" bookkeeping
    % variable (see pd_rule.m / pd_gmres.m: mCurrent tracks the growth
    % floor pd_rule falls back on when the PD law proposes m < mMin).
    m = mInitial;
    mCurrent = mInitial;

    % Sliding window of LGMRES-style error approximation vectors, newest
    % vector last.  nZ tracks how many columns are actually populated.
    zMat = zeros(n, l);
    nZ = 0;

    tic(); % start wall-clock timer

    % -------------------------------------------------------------------
    % Cycle 1: plain GMRES(m).  Neither augmentation strategy, nor the
    % restart-parameter growth law (which needs cycle history), applies
    % yet.
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

    zMat(:, 1) = zCycle;
    nZ = 1;

    stagnating = (relresvec(nCycles + 1) / relresvec(nCycles) >= normY);
    if stagnating
        % Cycle 1's basis is a plain (non-augmented) Arnoldi basis, so the
        % cheap H'-based formula for Fold is exact here -- see slgmres_e.m
        % for why this shortcut is only valid on this first cycle.
        Fold = H(1:s, 1:s)';
        G = Rs' * Rs;
        dy = harmonic_ritz_vectors(Fold, G, d, V, eigstol);
    end

    % =========================================================================
    % ----> Main loop: cycles 2, 3, ...
    % =========================================================================

    while flag == 0 && nCycles < maxit

        % ------------------------------------------------------------------
        % Control block: on a stagnating cycle, grow m via the PD law
        % (pd_rule, Algorithm 1 of [2]) before building this cycle's
        % subspace.  m is left untouched on non-stagnating (LGMRES-style)
        % cycles -- growth only happens in direct response to detected
        % slowdown, matching [1]'s Adaptive_PD_lgmres_e.m reference.
        % ------------------------------------------------------------------
        if stagnating
            miter = pd_rule(m, n, mCurrent, mMin, mMax, mStep, ...
                            relresvec, nCycles + 1, alphaP, alphaD);
            m = miter(1);
            mCurrent = miter(2);
        end

        r = b - A * x;
        beta = norm(r);
        v1 = r / beta;

        if ~stagnating
            % --- LGMRES-style cycle ---
            % See slgmres_e.m for why the fliplr placement here is the
            % opposite of the GMRES-E-style branch below.
            lUse = min(nZ, l);
            [H, V, s] = ...
                augmented_gram_schmidt_arnoldi(A, v1, m, zMat(:, 1:lUse));
            [HUpTri, g] = plane_rotations(H, beta);
            Rs = HUpTri(1:s, 1:s);
            gs = g(1:s);
            minimizer = Rs \ gs;
            V(:, m + 1:s) = fliplr(zMat(:, 1:lUse));
        else
            % --- GMRES-E-style cycle (using the just-grown m) ---
            %
            % dy IS sliced to d columns here, matching slgmres_e.m: this
            % family of algorithms is built on Cabral's reference
            % implementations (jcc_codigos_may_2023/), whose GMRES-E-style
            % branch (e.g. Adaptive_PD_lgmres_e.m) hard-codes s = m + d
            % and only ever reads dy(:, 1:d) -- functionally identical to
            % this explicit slice, even though harmonic_ritz_vectors can
            % return more columns when a harmonic Ritz value is complex.
            % See slgmres_e.m's GMRES-E-style branch for the full
            % rationale, and gmres_e.m (which has no such reference to
            % match) for why that solver instead uses dy in full.
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
        % Decide the augmentation strategy for the NEXT cycle -- see
        % slgmres_e.m for the full rationale, which carries over unchanged
        % here (the growth of m does not affect this decision).
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
