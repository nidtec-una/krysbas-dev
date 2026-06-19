function [x, flag, relresvec, kdvec, time] = ...
    gmres_dr(A, b, m, k, tol, maxit, xInitial, varargin)
    % GMRES-DR algorithm
    %
    %   Description:
    %   ------------
    %
    %   GMRES-DR ("GMRES with Deflated Restarting") is a restarted GMRES
    %   variant that deflates a k-dimensional invariant subspace at every
    %   restart by recycling the k harmonic Ritz vectors from the previous
    %   cycle directly into the starting basis of the next cycle.
    %
    %   Unlike GMRES-E(m, d), which augments a fresh Krylov subspace of
    %   dimension m with d extra eigenvectors (total dimension m + d),
    %   GMRES-DR(m, k) keeps the subspace dimension fixed at m per cycle:
    %   the first k basis vectors are the recycled Ritz vectors and the
    %   remaining m - k vectors come from a standard Arnoldi process started
    %   at the normalized, deflated residual.
    %
    %   A key efficiency property: the first k columns of the Hessenberg
    %   matrix for each restart cycle can be formed without any new
    %   matrix-vector products, reusing the Arnoldi factorization from the
    %   previous cycle (see subsection 2.3 of [1]).
    %
    %   Convergence is accelerated because the recycled Ritz vectors
    %   approximate the invariant subspace associated with the k smallest
    %   eigenvalues; deflating these eigenvalues reduces the effective
    %   condition number seen by the remaining GMRES iterations.
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
    %               Number of harmonic Ritz vectors to deflate and recycle
    %               at each restart.  Must satisfy 0 < k < m.
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
    %               1 if the relative residual satisfies relresvec(end) < tol
    %               within maxit cycles; 0 otherwise.
    %
    %   relresvec:  (cycles+1)-by-1 vector
    %               Relative residual norm after each cycle, starting from 1
    %               (the initial relative residual).
    %
    %   kdvec:      (cycles+1)-by-1 vector
    %               Krylov subspace dimension used at the corresponding cycle.
    %               For GMRES-DR this is m at every cycle (except possibly the
    %               last if happy breakdown occurs).
    %
    %   time:       float
    %               Wall-clock time in seconds.
    %
    %
    %   Algorithm Outline:
    %   ------------------
    %
    %   Cycle 1 (no Ritz vectors available yet):
    %     1. Arnoldi: build an (m+1)-by-m Hessenberg H and an n-by-m
    %        orthonormal basis V via modified_gram_schmidt_arnoldi.
    %     2. Least-squares: solve min ||beta*e1 - H*y|| via plane_rotations.
    %     3. Update: x = x0 + V * y.
    %     4. Convergence check.
    %     5. Ritz: compute k harmonic Ritz pairs via harmonic_ritz_vectors.
    %        Keep the Ritz vectors Yk = V * Ek (n-by-k) and the
    %        Arnoldi-space eigenvectors Ek (m-by-k).
    %
    %   Cycles 2, 3, ... (deflated restarts):
    %     6. QR: orthonormalize Yk -> [Vk, Rk] (thin QR).
    %        Vk (n-by-k) is the new deflation basis; Rk (k-by-k) is the
    %        upper-triangular factor.
    %     7. Residual: compute rc = b - A*x.  Orthogonalize rc against Vk
    %        (modified Gram-Schmidt) to obtain the deflated residual r_tilde
    %        and v_{k+1} = r_tilde / ||r_tilde||.
    %     8. Free columns of H (no new matrix-vector products):
    %        Using the previous Arnoldi factorization A*Vprev = Vprev_ext*Hprev,
    %        express A*Vk as a linear combination of Vprev_ext.  Then project
    %        onto the new basis [Vk, v_{k+1}] to fill H(:, 1:k).
    %        Details: let F = Ek / Rk  (m-by-k, coefficients of Vk in Vprev),
    %                     T = Hprev * Ek / Rk  ((m+1)-by-k),
    %           H(1:k,   1:k) = F' * T(1:m, :),
    %           H(k+1,   1:k) = g_tilde' * T / ||g_tilde||
    %        where g_tilde is the Arnoldi-space representation of r_tilde:
    %           g_res   = beta_prev*e1 - Hprev*yprev   (previous GMRES residual)
    %           g_tilde = [(I - F*F') * g_res(1:m);  g_res(m+1)]
    %     9. Arnoldi extension: run m - k standard Arnoldi steps from v_{k+1},
    %        orthogonalising against all of [Vk, v_{k+1}, ...] to fill
    %        H(:, k+1:m) and V(:, k+2:m+1).
    %    10. Least-squares: the right-hand side phi for this cycle is not
    %        beta*e1 but rather the projection of rc onto the new basis:
    %           phi(1:k)   = Vk' * rc  (from step 7's Gram-Schmidt coefficients)
    %           phi(k+1)   = ||r_tilde||
    %           phi(k+2:m+1) = 0
    %        Solve min ||phi - H*y|| using MATLAB's backslash on the full
    %        (m+1)-by-m system.
    %    11. Update: x = x + V * y.  Convergence check.
    %    12. Ritz: compute k new harmonic Ritz pairs for the next cycle.
    %
    %
    %   References:
    %   -----------
    %
    %   [1] Morgan, R. B. (2002). GMRES with deflated restarting.
    %   SIAM Journal on Scientific Computing, 24(1), 20-37.
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

    % =========================================================================
    % ----> Sanity checks on vector b
    % =========================================================================

    % Check whether b is non-empty
    if isempty(b)
        error("Vector b cannot be empty.");
    end

    % Check whether b is a column vector
    [rowsb, colsb] = size(b);
    if colsb ~= 1
        error("Vector b must be a column vector.");
    end

    % Check whether b has the correct number of rows
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

    % m == n: fall back to built-in unrestarted GMRES
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
    % ----> Default value and sanity checks for k
    % =========================================================================

    if (nargin < 4) || isempty(k)
        k = min(m - 1, 3);
    end

    % k == 0: no deflation, fall back to built-in restarted GMRES(m)
    if k == 0
        tic();
        [gmres_x, gmres_flag, ~, ~, resvec] = gmres(A, b, m);
        time = toc();
        x = gmres_x;
        flag = (gmres_flag == 0);
        relresvec = resvec ./ resvec(1);
        kdvec = m .* ones(length(relresvec), 1);
        return
    end

    % k must be strictly less than m so that at least one Arnoldi step is
    % performed per cycle beyond the deflation vectors.
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

    % Initial residual and its norm
    r0 = b - A * xInitial;
    res0 = norm(r0);

    % Relative residual history (entry 1 = before any iteration = 1.0)
    relresvec = zeros(maxit + 1, 1);
    relresvec(1) = 1.0;

    % Krylov dimension history (one entry per cycle)
    kdvec = zeros(maxit + 1, 1);

    nCycles = 0;  % number of completed cycles

    tic();  % start wall-clock timer

    % =========================================================================
    % ----> Cycle 1: standard GMRES(m) (no Ritz vectors available yet)
    % =========================================================================

    beta = res0;
    v1 = r0 / beta;

    % Arnoldi factorisation: A * V(:,1:s) = V(:,1:s+1) * H(1:s+1,1:s)
    % s may be less than m if happy breakdown occurs.
    [H, V, s] = modified_gram_schmidt_arnoldi(A, v1, m);

    % QR factorisation of H via Givens rotations; the rotated RHS g has
    % |g(s+1)| = ||r|| (the residual norm after this cycle).
    [HUpTri, g] = plane_rotations(H, beta);

    % Solve the triangular least-squares system
    Rs = HUpTri(1:s, 1:s);
    yprev = Rs \ g(1:s);

    % Update solution
    x = xInitial + V * yprev;

    nCycles = nCycles + 1;
    relresvec(nCycles + 1) = abs(g(s + 1)) / res0;
    kdvec(nCycles + 1) = s;

    % Check convergence after cycle 1
    if relresvec(nCycles + 1) < tol
        flag = 1;
        relresvec = relresvec(1:nCycles + 1);
        kdvec = kdvec(1:nCycles + 1);
        time = toc();
        return
    end

    % Store data from cycle 1 needed by subsequent deflated restarts:
    %   Hprev   : the (s+1)-by-s Hessenberg from the previous cycle
    %   Vprev   : the n-by-(s+1) Arnoldi basis (includes v_{s+1})
    %   g_res   : the GMRES residual in Arnoldi-space coordinates,
    %             defined as  g_res = beta*e1 - H*yprev
    %             Note: by construction, Vprev_ext * g_res = rc
    %   beta_eff: ||r0|| for cycle 1 (used to form g_res)
    Hprev = H;  % (s+1)-by-s, possibly truncated by happy breakdown
    % Vprev_ext = [V(:,1:s), v_{s+1}] where v_{s+1} is stored in
    % the (s+1)-th column only when s < m (happy breakdown gives V with s cols).
    % In the non-breakdown case modified_gram_schmidt_arnoldi returns V with s
    % columns (= m columns); the (s+1)-th Arnoldi vector is NOT stored in V
    % because it is not needed for the GMRES update.  We therefore have to
    % reconstruct it when s == m (no breakdown):
    %
    %   v_{m+1} = (A*v_m - V*H(:,m)) / H(m+1,m)
    %
    % This costs one extra matrix-vector product (1 matvec) at the end of
    % cycle 1 only.  All subsequent cycles use the efficient no-matvec
    % formula for the first k columns of H.
    if s == m
        % No happy breakdown: reconstruct v_{m+1} for the deflated restart.
        % A*v_m is already "inside" the Arnoldi recurrence but not stored.
        w_last = A * V(:, s);
        for i = 1:s
            w_last = w_last - H(i, s) * V(:, i);
        end
        % H(s+1, s) is the subdiagonal entry; v_{s+1} = w_last / H(s+1, s).
        v_ext = w_last / H(s + 1, s);
        Vprev_ext = [V, v_ext];  % n-by-(s+1)
    else
        % Happy breakdown: V already contains s columns; H(s+1,s) == 0 and
        % there is no (s+1)-th Arnoldi vector (the Krylov subspace is
        % invariant).  In this case, the algorithm should have converged.
        Vprev_ext = V;  % n-by-s (exceptional, convergence was missed above)
    end

    % Arnoldi-space residual from cycle 1:
    %   g_res is an (s+1)-by-1 vector satisfying  Vprev_ext * g_res = rc.
    g_res = zeros(s + 1, 1);
    g_res(1) = beta;
    g_res = g_res - Hprev * yprev;  % = beta*e1 - H*y

    % Harmonic Ritz computation for the first deflated restart.
    % F_eig  : m-by-k matrix whose columns are the eigenvectors of the
    %          k-smallest harmonic Ritz problem in Arnoldi-space coordinates.
    % The existing helper harmonic_ritz_vectors(F, G, k, V, tol) needs:
    %   F = H(1:s,1:s)'   (the square part transposed)
    %   G = Rs'*Rs         (Gram matrix from the least-squares QR)
    %   k = number of Ritz vectors wanted
    %   V = the Arnoldi basis (for forming n-space Ritz vectors)
    Fsq = Hprev(1:s, 1:s)';
    G = Rs' * Rs;
    % harmonic_ritz_vectors returns dy: n-by-k matrix of Ritz vectors in
    % the physical (n-dimensional) space, and implicitly uses Ek internally.
    % We need both the n-space vectors AND the Arnoldi-space eigenvectors Ek
    % to build the deflated Hessenberg efficiently.  Because the existing
    % function does not return Ek, we re-solve the generalized eigenvalue
    % problem here to obtain it.
    %
    % Generalized eigenvalue problem (harmonic Ritz, [1] p. 1161, step 5):
    %   Fsq * u = lambda * G * u
    % Take the k eigenpairs with smallest |lambda|; columns of Ek are the
    % eigenvectors in Arnoldi space; n-space Ritz vectors are Vk = V * Ek.
    opts_eig.tol = tol;
    opts_eig.v0 = ones(s, 1);
    [Ek_raw, Dk] = eigs(Fsq, G, k, 'LM', opts_eig);

    % Sort by ascending |eigenvalue| to get the k smallest
    [~, idx] = sort(abs(diag(Dk)));
    % Take real part: for a real problem the imaginary components of the
    % harmonic Ritz vectors are numerical noise introduced by eigs.
    Ek = real(Ek_raw(:, idx));  % s-by-k, Arnoldi-space eigenvectors

    % n-space Ritz vectors (before orthonormalization)
    Yk = V * Ek;  % n-by-k

    % =========================================================================
    % ----> Deflated restart cycles (cycles 2, 3, ...)
    % =========================================================================

    for cycle = 2:maxit

        % ------------------------------------------------------------------
        % Step 1: QR-orthonormalize the k Ritz vectors.
        %
        % Vk (n-by-k, orthonormal columns) will serve as the first k basis
        % vectors of this cycle's subspace.  Rk_qr (k-by-k upper triangular)
        % satisfies  Yk = Vk * Rk_qr.  Consequently the Arnoldi-space
        % representation of Vk in terms of the PREVIOUS Arnoldi basis Vprev is:
        %
        %   Vk = Vprev(:,1:s) * Ek_norm
        %
        % where  Ek_norm = Ek / Rk_qr  (s-by-k).
        % ------------------------------------------------------------------
        [Vk, Rk_qr] = qr(Yk, 0);  % thin QR: Vk is n-by-k, Rk_qr is k-by-k

        % Coefficient matrix relating Vk to the previous Arnoldi basis
        % (used below to form H's first k columns without new matvecs)
        Ek_norm = Ek / Rk_qr;  % s-by-k

        % ------------------------------------------------------------------
        % Step 2: Compute the current residual and deflate it.
        %
        % rc = b - A*x  (one matrix-vector product)
        %
        % Orthogonalize rc against Vk using modified Gram-Schmidt to obtain
        % the deflated residual r_tilde and the (k+1)-th basis vector v_{k+1}.
        %
        % The Gram-Schmidt coefficients  hgs(i) = Vk(:,i)' * rc  are the
        % first k entries of the RHS vector phi used in the least-squares
        % problem (step 10 of the algorithm outline above).
        % ------------------------------------------------------------------
        rc = b - A * x;

        % Modified Gram-Schmidt projection
        r_tilde = rc;
        hgs = zeros(k, 1);  % Gram-Schmidt coefficients (= phi(1:k))
        for i = 1:k
            hgs(i) = Vk(:, i)' * r_tilde;
            r_tilde = r_tilde - hgs(i) * Vk(:, i);
        end

        norm_r_tilde = norm(r_tilde);

        % Check for happy breakdown: rc is already in span(Vk)
        if norm_r_tilde < eps * norm(rc)
            % The residual lies (numerically) in the deflation subspace.
            % The current x is the best approximation we can achieve with
            % the deflation subspace alone; report relative residual and exit.
            relresvec(nCycles + 2) = norm(rc) / res0;
            kdvec(nCycles + 2) = k;
            nCycles = nCycles + 1;
            flag = (relresvec(nCycles + 1) < tol);
            relresvec = relresvec(1:nCycles + 1);
            kdvec = kdvec(1:nCycles + 1);
            time = toc();
            return
        end

        v_kp1 = r_tilde / norm_r_tilde;  % v_{k+1}: (k+1)-th basis vector

        % ------------------------------------------------------------------
        % Step 3: Form the first k columns of the new Hessenberg matrix H
        %         WITHOUT any new matrix-vector products.
        %
        % From the previous cycle's Arnoldi factorisation:
        %   A * Vprev(:,1:s) = Vprev_ext * Hprev       (*)
        % where Vprev_ext = [Vprev(:,1:s), v_{s+1}]  (n-by-(s+1)).
        %
        % The new deflation basis Vk satisfies:
        %   Vk = Vprev(:,1:s) * Ek_norm
        % so substituting into (*):
        %   A * Vk = Vprev_ext * (Hprev * Ek_norm)    (**)
        %
        % Define  T = Hprev * Ek_norm   ((s+1)-by-k, cheap to compute).
        % Then  A*Vk = Vprev_ext * T.
        %
        % The (i,j) entry of H(:, 1:k) is  v_i' * (A*v_j):
        %
        % For i = 1,...,k (deflation vectors):
        %   v_i = Vprev(:,1:s) * Ek_norm(:,i)
        %   v_i' * Vprev_ext = Ek_norm(:,i)' * [I_s, 0]  (only first s entries)
        %   => H(i, j) = Ek_norm(:,i)' * T(1:s, j) = (Ek_norm' * T(1:s,:))(i,j)
        %
        % For i = k+1 (the deflated residual vector v_{k+1}):
        %   v_{k+1}' * Vprev_ext = ?
        %   We know  rc = Vprev_ext * g_res  (from step 7 of cycle 1, or the
        %   analogous formula maintained across cycles), so the Arnoldi-space
        %   image of rc under the previous basis is g_res.
        %   The orthogonalised deflated residual r_tilde in Arnoldi space is:
        %     g_tilde = [(I - Ek_norm*Ek_norm') * g_res(1:s);  g_res(s+1)]
        %   and  v_{k+1} = r_tilde / norm_r_tilde
        %                = Vprev_ext * g_tilde / ||g_tilde||.
        %   Therefore:
        %     v_{k+1}' * Vprev_ext = g_tilde' / norm_r_tilde
        %   => H(k+1, j) = g_tilde' * T(:, j) / norm_r_tilde
        % ------------------------------------------------------------------

        % T = Hprev * Ek_norm  ((s+1)-by-k)
        T = Hprev * Ek_norm;

        % Top k-by-k block: H(1:k, 1:k)
        H_def = zeros(m + 1, m);  % full Hessenberg for this cycle
        H_def(1:k, 1:k) = Ek_norm' * T(1:s, :);

        % Arnoldi-space deflated residual (s+1-vector):
        %   g_tilde(1:s)   = (I - Ek_norm*Ek_norm') * g_res(1:s)
        %   g_tilde(s+1)   = g_res(s+1)
        g_tilde = g_res;
        g_tilde(1:s) = g_res(1:s) - Ek_norm * (Ek_norm' * g_res(1:s));

        % Row k+1 of H's first k columns
        H_def(k + 1, 1:k) = (g_tilde' * T) / norm_r_tilde;

        % ------------------------------------------------------------------
        % Step 4: Assemble the full Arnoldi basis for this cycle.
        %
        % The basis V_cycle is n-by-(m+1):
        %   V_cycle(:, 1:k)     = Vk   (deflation vectors, already computed)
        %   V_cycle(:, k+1)     = v_{k+1}  (deflated residual direction)
        %   V_cycle(:, k+2:m+1) = filled by Arnoldi extension below
        % ------------------------------------------------------------------
        V_cycle = zeros(n, m + 1);
        V_cycle(:, 1:k) = Vk;
        V_cycle(:, k + 1) = v_kp1;

        % ------------------------------------------------------------------
        % Step 5: Arnoldi extension from v_{k+1} for m - k additional steps.
        %
        % At step j (j = k+1, ..., m), we compute:
        %   w = A * V_cycle(:, j)
        %   Orthogonalise w against V_cycle(:, 1:j)   (all previous vectors)
        %   H(1:j+1, j) = Gram-Schmidt coefficients
        %   V_cycle(:, j+1) = w / H(j+1, j)
        %
        % This fills H_def(:, k+1:m) and V_cycle(:, k+2:m+1).
        % ------------------------------------------------------------------
        s_new = m;  % will be updated on happy breakdown
        for j = k + 1:m
            w = A * V_cycle(:, j);
            for i = 1:j
                H_def(i, j) = w' * V_cycle(:, i);
                w = w - H_def(i, j) * V_cycle(:, i);
            end
            H_def(j + 1, j) = norm(w);

            if H_def(j + 1, j) == 0
                % Happy breakdown: Krylov subspace is invariant.
                s_new = j;
                V_cycle = V_cycle(:, 1:s_new);
                H_def = H_def(1:s_new + 1, 1:s_new);
                break
            else
                V_cycle(:, j + 1) = w / H_def(j + 1, j);
            end
        end

        % ------------------------------------------------------------------
        % Step 6: Solve the least-squares problem for this cycle.
        %
        % The right-hand side phi is NOT beta*e1 (as in standard GMRES) but
        % instead the projection of rc onto the new basis [Vk, v_{k+1}, ...]:
        %
        %   phi(1:k)        = hgs      (Gram-Schmidt coefficients from step 2)
        %   phi(k+1)        = norm_r_tilde  (norm of deflated residual)
        %   phi(k+2:s_new+1) = 0        (orthogonality from Arnoldi)
        %
        % Solve:  min_{y} || phi - H_def(1:s_new+1, 1:s_new) * y ||
        %
        % We use MATLAB's backslash (\) on the rectangular system because the
        % RHS phi is not of the beta*e1 form required by plane_rotations.
        % ------------------------------------------------------------------
        phi = zeros(s_new + 1, 1);
        phi(1:k) = hgs;
        phi(k + 1) = norm_r_tilde;

        y = H_def(1:s_new + 1, 1:s_new) \ phi;

        % ------------------------------------------------------------------
        % Step 7: Update the approximate solution.
        % ------------------------------------------------------------------
        x = x + V_cycle(:, 1:s_new) * y;

        nCycles = nCycles + 1;

        % Relative residual for this cycle (exact residual norm)
        r_new = b - A * x;
        relresvec(nCycles + 1) = norm(r_new) / res0;
        kdvec(nCycles + 1) = s_new;

        % Check convergence
        if relresvec(nCycles + 1) < tol
            flag = 1;
            relresvec = relresvec(1:nCycles + 1);
            kdvec = kdvec(1:nCycles + 1);
            time = toc();
            return
        end

        % ------------------------------------------------------------------
        % Step 8: Compute k new harmonic Ritz pairs for the next cycle.
        %
        % The generalized eigenvalue problem is the same as in GMRES-E:
        %   H(1:s,1:s)' * u = lambda * (Rs_new' * Rs_new) * u
        % where Rs_new comes from a QR factorisation of H_def.
        %
        % Rather than calling the existing harmonic_ritz_vectors helper
        % (which does not return Ek), we repeat the eigenvalue computation
        % here to obtain both the n-space Ritz vectors (Yk) and the
        % basis-space eigenvectors (Ek) needed for the next cycle's free
        % columns of H.
        % ------------------------------------------------------------------

        % QR factorisation of the Hessenberg (for G = Rs'*Rs)
        [~, Rs_new] = qr(H_def(1:s_new + 1, 1:s_new), 0);
        Rs_new = Rs_new(1:s_new, :);

        Fsq_new = H_def(1:s_new, 1:s_new)';
        G_new = Rs_new' * Rs_new;

        opts_eig.v0 = ones(s_new, 1);
        [Ek_raw_new, Dk_new] = eigs(Fsq_new, G_new, k, 'LM', opts_eig);
        [~, idx_new] = sort(abs(diag(Dk_new)));
        % Take real part: imaginary components are numerical noise from eigs.
        Ek = real(Ek_raw_new(:, idx_new));  % s_new-by-k

        % n-space Ritz vectors for next cycle
        Yk = V_cycle(:, 1:s_new) * Ek;  % n-by-k

        % ------------------------------------------------------------------
        % Step 9: Store previous-cycle data needed for the next deflated
        %         restart's free-column computation.
        %
        % Hprev : Hessenberg from this cycle  ((s_new+1)-by-s_new)
        % g_res : Arnoldi-space image of the new residual r_new.
        %         Satisfies  V_cycle * g_res = b - A*x_new  (= r_new),
        %         which follows from the Arnoldi relation
        %           A*V_cycle(:,1:s_new) = V_cycle * H_def
        %         and from  r_new = V_cycle*phi - V_cycle*(H_def*y).
        %         V_cycle itself need NOT be stored: the free-column formula
        %         for the NEXT cycle uses only Hprev, Ek_norm, g_res, and s.
        % s     : subspace dimension of this cycle (= s_new, used as array
        %         size in the next cycle's Arnoldi-space computations)
        % ------------------------------------------------------------------
        Hprev = H_def;
        g_res = phi - H_def(1:s_new + 1, 1:s_new) * y;
        s = s_new;

    end  % deflated restart loop

    % =========================================================================
    % ----> Maximum number of cycles reached without convergence
    % =========================================================================

    relresvec = relresvec(1:nCycles + 1);
    kdvec = kdvec(1:nCycles + 1);
    time = toc();

end
