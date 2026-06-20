function [x, flag, relresvec, kdvec, time, stats] = ...
    gmres_dr(A, b, m, k, tol, maxit, xInitial, varargin)
    % GMRES-DR algorithm
    %
    %   Description:
    %   ------------
    %
    %   GMRES-DR ("GMRES with Deflated Restarting") is a restarted GMRES
    %   variant that recycles k harmonic Ritz vectors across restart cycles
    %   to accelerate convergence by deflating the k smallest eigenvalues.
    %
    %   Unlike GMRES-E(m, d), which augments a fresh m-step Krylov subspace
    %   with d extra eigenvectors (total dimension m + d), GMRES-DR(m, k)
    %   keeps the subspace dimension fixed at m per cycle: the first keep
    %   basis vectors are the recycled Schur vectors (keep <= k, see below)
    %   and the remaining m - keep vectors come from standard Arnoldi.
    %
    %   This implementation follows Morgan (2002) and uses:
    %
    %   (a) INCREMENTAL QR (qrupdate_gs): the thin QR of the growing
    %       Hessenberg H is updated column by column via modified
    %       Gram-Schmidt.  This avoids refactorising H at every step and
    %       provides a cheap residual-norm estimate at each inner iteration.
    %
    %   (b) SCHUR-BASED HARMONIC RITZ: the harmonic Ritz problem is solved
    %       via the explicit harmonic matrix
    %
    %         Ht = H_sq + h_{m+1,m}^2 * (H_sq' \ e_m) * e_m'
    %
    %       followed by a real Schur decomposition and ordschur reordering.
    %       This is numerically more stable than the generalized eigenvalue
    %       formulation eigs(F, G, k, 'LM') used in GMRES-E, which can fail
    %       when G = Rs'*Rs is ill-conditioned.
    %
    %   (c) THICK RESTART via orthogonal basis rotation (Morgan 2002,
    %       Wu & Simon 1999): the k Schur vectors Pk and the residual
    %       direction rc are combined into a new (keep+1)-column basis
    %       Pkp1 via thin QR.  The Hessenberg and V are then updated as
    %
    %         H  <-- Pkp1' * H * Pk     (now (keep+1)-by-keep)
    %         V  <-- V * Pkp1           (rotate basis in-place)
    %         Vr <-- Pkp1' * rc         (residual in new coordinates)
    %
    %       This is a pure orthogonal transformation: no matrix-vector
    %       products with A are needed for the recycle step, and no
    %       approximation of A*Vk from a previous factorisation is used.
    %
    %   (d) 2x2 SCHUR BLOCK HANDLING: when the (k+1,k) entry of the
    %       ordered Schur form T is non-zero, a complex conjugate pair
    %       straddles the cut.  We set keep = k+1 to avoid splitting the
    %       pair, and reduce back to k at the next restart.
    %
    %   Signature:
    %   ----------
    %
    %   [x, flag, relresvec, kdvec, time] = ...
    %       gmres_dr(A, b, m, k, tol, maxit, xInitial)
    %
    %
    %   Input Parameters:
    %   -----------------
    %
    %   A:          n-by-n matrix (dense or sparse)
    %               Left-hand side of the linear system Ax = b.
    %
    %   b:          n-by-1 vector
    %               Right-hand side of the linear system Ax = b.
    %
    %   m:          int, optional
    %               Subspace dimension per restart cycle.  If m == n,
    %               the built-in unrestarted GMRES is used.  Must satisfy
    %               1 <= m <= n.  Default: min(n, 10).
    %
    %   k:          int, optional
    %               Number of harmonic Ritz vectors to recycle at each
    %               restart.  Must satisfy 0 < k < m.
    %               Default: min(m - 1, 3).
    %               Setting k = 0 falls back to standard restarted GMRES(m).
    %
    %   tol:        float, optional
    %               Relative residual tolerance.  The algorithm stops when
    %               ||r||/||r0|| < tol.  Default: 1e-6.
    %
    %   maxit:      int, optional
    %               Maximum number of restart cycles.  Default: min(n, 10).
    %
    %   xInitial:   n-by-1 vector, optional
    %               Initial guess.  Default: zeros(n, 1).
    %
    %
    %   Output Parameters:
    %   ------------------
    %
    %   x:          n-by-1 vector
    %               Approximate solution to Ax = b.
    %
    %   flag:       integer (0 or 1)
    %               1 if ||r||/||r0|| < tol within maxit cycles; 0 otherwise.
    %
    %   relresvec:  (cycles+1)-by-1 vector
    %               Relative residual norm after each cycle, starting from 1.
    %               Estimated from the Arnoldi-space residual norm, which
    %               equals the true residual norm when the Arnoldi relation
    %               A*V = V_ext*H holds exactly.
    %
    %   kdvec:      (cycles+1)-by-1 vector
    %               Krylov subspace dimension used at each cycle.  Equals m
    %               at every full cycle; may be smaller at early exit.
    %
    %   time:       float
    %               Wall-clock time in seconds.
    %
    %   stats:      struct (optional, only populated when requested)
    %               Diagnostic counters for the harmonic Ritz path chosen
    %               at each restart cycle:
    %                 stats.n_eigs  : number of cycles that used eigs
    %                                 (well-conditioned G, rcond > 1e-14)
    %                 stats.n_schur : number of cycles that used the Schur
    %                                 fallback (ill-conditioned G)
    %
    %
    %   References:
    %   -----------
    %
    %   [1] Morgan, R. B. (2002). GMRES with deflated restarting.
    %   SIAM Journal on Scientific Computing, 24(1), 20-37.
    %
    %   [2] Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for
    %   large symmetric eigenvalue problems. SIAM Journal on Matrix
    %   Analysis and Applications, 22(2), 602-616.
    %
    %   [3] GMRES-SDR reference implementation by the NLA Group, University
    %   of Manchester.  The Schur-based thick restart and incremental QR
    %   update adopted in this file are adapted from the gmres_dr.m and
    %   qrupdate_gs.m files in that repository.
    %   https://github.com/nla-group/GMRES-SDR
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
    elseif nargin > 7
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
        stats = struct('n_eigs', 0, 'n_schur', 0);
        return
    end

    if m > n
        error("m must satisfy: 1 <= m <= n.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for k
    % =========================================================================

    if (nargin < 4) || isempty(k)
        k = min(m - 1, 3);
    end

    if k == 0
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b, m);
        time = toc();
        x = gmres_x;
        flag = (gmres_flag == 0);
        relresvec = resvec ./ resvec(1);
        kdvec = m .* ones(length(relresvec), 1);
        stats = struct('n_eigs', 0, 'n_schur', 0);
        return
    end

    if k >= m
        error("k must satisfy: 0 < k < m.");
    end

    % =========================================================================
    % ----> Default value and sanity checks for tol
    % =========================================================================

    if (nargin < 5) || isempty(tol)
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

    if (nargin < 6) || isempty(maxit)
        maxit = min(n, 10);
    end

    % =========================================================================
    % ----> Default value and sanity checks for xInitial
    % =========================================================================

    if (nargin < 7) || isempty(xInitial)
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
    % ----> GMRES-DR Algorithm starts here
    % =========================================================================

    flag = 0;
    x    = xInitial;
    r0   = b - A * x;
    beta = norm(r0);    % initial residual norm (used for relative scaling)

    % Relative residual and dimension history (one entry per cycle)
    relresvec = zeros(maxit + 1, 1);
    relresvec(1) = 1.0;
    kdvec   = zeros(maxit + 1, 1);
    n_cycles       = 0;
    n_eigs_cycles  = 0;   % cycles that used the eigs primary path
    n_schur_cycles = 0;   % cycles that used the Schur fallback

    % ------------------------------------------------------------------
    % Arnoldi basis V: pre-allocated for m + 1 columns.  After each
    % thick restart, columns 1:keep+1 are overwritten; the rest are
    % filled by the next Arnoldi extension.
    % Vr: coordinates of the current residual in the V-basis.
    %     Initially just [beta] (single scalar, since V(:,1) = r/beta).
    % H:  Hessenberg matrix, grown column-by-column during Arnoldi and
    %     shrunk to (keep+1)-by-keep after each thick restart.
    % keep: number of Schur vectors recycled from the previous cycle.
    %       Starts at 0 (no recycling on the first cycle).
    % ------------------------------------------------------------------
    V  = zeros(n, m + 1);
    V(:, 1) = r0 / beta;
    vr   = beta;        % (j+1)-vector grown as Arnoldi proceeds
    H    = zeros(0, 0); % empty; grows to (m+1)-by-m over the cycle
    keep = 0;           % no recycled vectors yet

    tic();  % start wall-clock timer

    for cycle = 1:maxit

        % ------------------------------------------------------------------
        % Each cycle targets a total subspace of dimension m.
        %   - Columns 1:keep of V are the recycled Schur vectors (already
        %     set by the previous thick restart, or empty on cycle 1).
        %   - Columns keep+1:m are added by the Arnoldi extension below.
        %   - Column m+1 is a temporary workspace for the last Arnoldi
        %     vector (not used in the least-squares solve).
        % ------------------------------------------------------------------

        % Initialise the incremental QR for the carry-over H block.
        % On cycle 1 (keep=0), q_qr and r_qr start empty.
        % On later cycles, H(1:keep+1, 1:keep) holds the rotated block
        % from the previous thick restart; qrupdate_gs processes it first.
        q_qr = zeros(keep + 1, 0);
        r_qr = zeros(0, 0);
        if keep > 0
            [q_qr, r_qr] = qrupdate_gs(H, q_qr, r_qr);
        end

        % Allocate (m+1)-by-m Hessenberg for this cycle, copying the
        % carry-over block into the top-left corner.
        H_cycle = zeros(m + 1, m);
        H_cycle(1:size(H, 1), 1:size(H, 2)) = H;
        H = H_cycle;

        % ------------------------------------------------------------------
        % Arnoldi extension: add m - keep new basis vectors.
        %
        % At each step j the loop:
        %   1. Applies A to obtain a new candidate w = A * V(:, j).
        %   2. Orthogonalises w against all previous V(:, 1:j) via MGS,
        %      filling column j of H.
        %   3. Updates the thin QR of H(1:j+1, 1:j) via qrupdate_gs.
        %   4. Solves the incremental least-squares problem
        %        min_{d} ||vr_c - H(1:j+1,1:j)*d||
        %      where vr_c is vr extended to length j+1 (zero-padded).
        %   5. Estimates the residual norm as ||rc|| = ||vr_c - H*d||.
        %      This equals the true residual norm because the Arnoldi
        %      relation gives ||b - A*(x+V*d)|| = ||V_ext * rc|| = ||rc||.
        %   6. Exits early if convergence is detected.
        % ------------------------------------------------------------------
        d  = zeros(m, 1);   % least-squares solution (updated each step)
        rc = vr;            % Arnoldi-space residual (updated each step)
        res = norm(vr);     % current residual norm estimate

        for j = keep + 1:m
            w = A * V(:, j);

            % Modified Gram-Schmidt orthogonalisation against all V(:,1:j)
            for i = 1:j
                H(i, j) = V(:, i)' * w;
                w = w - H(i, j) * V(:, i);
            end
            H(j + 1, j) = norm(w);
            V(:, j + 1)  = w / H(j + 1, j);

            % Extend vr with a zero so it matches the new (j+1) dimension.
            % This reflects that V(:,j+1) is orthogonal to all of r.
            vr_c = [vr; zeros(j + 1 - length(vr), 1)]; %#ok<AGROW>

            % Incremental QR update: add column j to the factorisation
            [q_qr, r_qr] = qrupdate_gs(H(1:j + 1, 1:j), q_qr, r_qr);

            % Least-squares solve via the updated QR
            d  = r_qr \ (q_qr' * vr_c);
            rc = vr_c - H(1:j + 1, 1:j) * d;
            res = norm(rc);

            % Convergence check (relative to initial residual beta)
            if res < tol * beta
                x = x + V(:, 1:j) * d;
                n_cycles = n_cycles + 1;
                relresvec(n_cycles + 1) = res / beta;
                kdvec(n_cycles + 1)     = j;
                flag = 1;
                relresvec = relresvec(1:n_cycles + 1);
                kdvec     = kdvec(1:n_cycles + 1);
                stats = struct('n_eigs', n_eigs_cycles, ...
                               'n_schur', n_schur_cycles);
                time = toc();
                return
            end
        end

        % ------------------------------------------------------------------
        % End of Arnoldi loop: update the solution and record history.
        % ------------------------------------------------------------------
        x = x + V(:, 1:m) * d;
        n_cycles = n_cycles + 1;
        relresvec(n_cycles + 1) = res / beta;
        kdvec(n_cycles + 1)     = m;

        if relresvec(n_cycles + 1) < tol
            flag = 1;
            relresvec = relresvec(1:n_cycles + 1);
            kdvec     = kdvec(1:n_cycles + 1);
            stats = struct('n_eigs', n_eigs_cycles, ...
                           'n_schur', n_schur_cycles);
            time = toc();
            return
        end

        % ------------------------------------------------------------------
        % Harmonic Ritz computation -- hybrid eigs / Schur strategy.
        %
        % G = R'*R  (from the incremental QR already computed this cycle)
        % F = H_sq' (square part of the Hessenberg, transposed)
        %
        % PRIMARY PATH -- eigs(F, G, k, 'LM'):
        %   Identical to the GMRES-E approach.  Fast and accurate when G is
        %   well-conditioned.  We sort the returned eigenvalues in ascending
        %   |lambda| order and QR-orthonormalise the eigenvectors to form Pk.
        %   real() is applied to strip conjugate-noise artefacts.
        %
        % FALLBACK PATH -- Schur decomposition of the harmonic matrix Ht:
        %   Used when rcond(G) < 1e-14 (suggested by J.C. Cabral), i.e. when
        %   G is too ill-conditioned for the generalized eigensolver.  We form
        %
        %     Ht = H_sq + h_{m+1,m}^2 * (H_sq' \ e_m) * e_m'
        %
        %   and compute the real Schur decomposition.  ordschur reorders the
        %   k Schur vectors with smallest |eigenvalue| to the front.  The
        %   Schur form naturally handles complex conjugate pairs as 2x2 blocks
        %   in T: if T(k+1,k) ~= 0 the cut falls inside such a block, so we
        %   set keep = k+1 to avoid splitting the pair.
        % ------------------------------------------------------------------
        g_mat = r_qr' * r_qr;           % G from the already-computed QR
        f_mat = H(1:m, 1:m)';           % F = H_sq'

        if rcond(g_mat) > 1e-14
            % --- eigs path (primary) ---
            n_eigs_cycles = n_eigs_cycles + 1;
            opts_eig.tol = tol;
            opts_eig.v0  = ones(m, 1);
            [ek_raw, dk] = eigs(f_mat, g_mat, k, 'LM', opts_eig);
            [~, idx] = sort(abs(diag(dk)));
            ek = real(ek_raw(:, idx));
            [pk, ~] = qr(ek, 0);        % orthonormalise eigenvectors
            keep = k;
        else
            % --- Schur path (fallback when G is ill-conditioned) ---
            n_schur_cycles = n_schur_cycles + 1;
            h_sub = H(m + 1, m);
            emk   = zeros(m, 1);
            emk(m) = 1;
            ht = f_mat' + h_sub^2 * (f_mat \ emk) * emk';
            [s_vecs, s_vals] = schur(ht);
            ritz = ordeig(s_vals);
            [~, ind] = sort(abs(ritz), 'ascend');
            sel = false(m, 1);
            sel(ind(1:k)) = true;
            [s_vecs, s_vals] = ordschur(s_vecs, s_vals, sel);
            if k > 0 && s_vals(k + 1, k) ~= 0
                keep = k + 1;
            else
                keep = k;
            end
            pk = s_vecs(:, 1:keep);
        end

        % ------------------------------------------------------------------
        % Thick restart: orthogonal basis rotation.
        %
        % pk (m-by-keep) is set above: either the QR-orthonormalised eigs
        % eigenvectors (primary path) or the ordered Schur vectors (fallback).
        % We construct the new (m+1)-by-(keep+1) basis Pkp1 by:
        %   1. Extending Pk with a zero bottom row: [(m-by-keep); (1-by-keep)]
        %   2. Appending rc (the current Arnoldi-space residual vector):
        %      the last column of Pkp1 will be the normalized residual
        %      direction after removing its Pk components.
        %   3. Thin QR of the result to make Pkp1 orthonormal.
        %   4. SIGN CORRECTION: restore the first keep columns to Pk exactly.
        %      The QR can flip column signs arbitrarily; without this step
        %      the recycled H block would have sign errors that corrupt
        %      convergence.
        %
        % After this, the thick-restart updates are:
        %   H  <-- Pkp1' * H(1:m+1, 1:m) * Pk    (now (keep+1)-by-keep)
        %   V(:,1:keep+1) <-- V(:,1:m+1) * Pkp1   (rotate basis)
        %   vr <-- Pkp1' * rc                      (residual coordinates)
        %
        % The V update overwrites only columns 1:keep+1; columns keep+2:m+1
        % are filled by the Arnoldi extension in the next cycle.
        % ------------------------------------------------------------------
        pkp1 = [pk; zeros(1, keep)];     % extend Pk with a zero bottom row
        pkp1 = [pkp1, rc];               % append residual direction
        [pkp1, ~] = qr(pkp1, 0);         % orthonormalize
        pkp1(1:m, 1:keep) = pk;          % restore sign of recycled columns

        H = pkp1' * H(1:m + 1, 1:m) * pk;
        V(:, 1:keep + 1) = V(:, 1:m + 1) * pkp1;
        vr = pkp1' * rc;

    end  % restart loop

    % =========================================================================
    % ----> Maximum number of cycles reached without convergence
    % =========================================================================

    relresvec = relresvec(1:n_cycles + 1);
    kdvec     = kdvec(1:n_cycles + 1);
    stats = struct('n_eigs', n_eigs_cycles, 'n_schur', n_schur_cycles);
    time = toc();

end
