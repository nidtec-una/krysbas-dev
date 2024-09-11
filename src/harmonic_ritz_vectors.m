function dy = harmonic_ritz_vectors(F, G, k, V, tol)
    % Harmonic Ritz Vectors function
    %
    %   This function is a modified implementation of the Rayleigh-Ritz method
    %   that finds good approximations to the smallest eigenvalues and its
    %   associated eigenvectors, and appends these ones to the next search
    %   subspaces. Refer to [1] for more information about this technique.
    %
    %   Signature:
    %   ----------
    %
    %   dy = harmonic_ritz_vectors(F, G, k, V, tol)
    %
    %
    %   Input Parameters:
    %   -----------------
    %
    %   F:          s-by-s maxtrix
    %               Matrix F from the generalized eigenvalue problem, as
    %               employed in step 5, p. 1161 of [1].
    %
    %   G:          s-by-1 vector
    %               Matrix G from the generalized eigenvalue problem, as
    %               employed in step 5, p. 1161 of [1].
    %
    %   k:          int
    %               Number of eigenvectors corresponding to a few of the
    %               smallest eigenvalues in magnitude. According to [1],
    %               "even just a few eigenvectors can make a big difference
    %               if the matrix has both small and large eigenvalues".
    %
    %   V:          n-by-s matrix
    %               A matrix from the modified Gram-Schmidt Arnoldi method,
    %               whose column vectors span the search subspace.
    %
    %   tol:        float
    %               Convergence tolerance for the built-in 'eigs' algorithm.
    %
    %
    %   Output parameters:
    %   ------------------
    %
    %   dy:         n-by-k matrix
    %               Matrix whose column vectors are the approximate
    %               eigenvectors of the matrix A.
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
    opts.tol = tol;
    s = size(F, 1);
    opts.v0 = ones(s, 1);
    E = zeros(s, k);
    D = zeros(k, 1);

    % Compute the approximate eigenvectors
    [E2, D2] = eigs(F, G, k, 'LM', opts);
    for p = 1:k
        D(p, 1) = abs(D2(p, p));
    end
    [~, I] = sort(D, 1);
    for q = 1:k
        E(:, q) = E2(:, I(q, 1));
    end
    dy0 = V * E; % Implements yi = Q * gi, step 5, p. 1161 of [1]

    dy = [];

    % If dy0 has complex components, its complex conjugate is also
    % an approximate eigenvector, hence we separate each eigenvector
    % into its real and complex parts and treat as two distinct vectors
    if isreal(dy0) == 0
        ij = 1;
        jj = 0;
        while size(dy, 2) <= k && ij <= k
            if isreal(dy0(:, ij)) == 0 && norm(real(dy0(:, ij))) > 0
                dy(:, jj + 1) = real(dy0(:, ij));
                jj = size(dy, 2);
                if ij <= k
                    dy(:, jj + 1) = abs(imag(dy0(:, ij)) * sqrt(1));
                    jj = size(dy, 2);
                    if ij < k
                        if dy0(:, ij) == conj(dy0(:, ij + 1))
                            ij = ij + 2;
                        else
                            ij = ij + 1;
                        end
                    end
                end
            else
                dy(:, jj + 1) = dy0(:, ij);
                ij = ij + 1;
                jj = size(dy, 2);
            end
        end
    else
        dy = dy0;
    end
end
