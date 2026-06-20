function [q_out, r_out] = qrupdate_gs(a, q_in, r_in)
    % Incremental modified Gram-Schmidt QR factorisation update.
    %
    %   Description:
    %   ------------
    %
    %   Given a tall matrix A and an existing partial thin QR factorisation
    %   of its first k columns, extends the factorisation to cover all n
    %   columns of A.  Each new column is orthogonalised against the running
    %   basis q_in using two passes of modified Gram-Schmidt (double
    %   reorthogonalisation) for numerical stability.
    %
    %   This utility is used by gmres_dr to maintain the QR of the growing
    %   Hessenberg matrix without recomputing it from scratch at every step.
    %
    %   Signature:
    %   ----------
    %
    %   [q_out, r_out] = qrupdate_gs(a, q_in, r_in)
    %
    %   Input Parameters:
    %   -----------------
    %
    %   a:      m-by-n matrix, m > n
    %           Full matrix whose thin QR is being built incrementally.
    %           Columns 1:k are assumed to be already factorised.
    %
    %   q_in:   m-by-k orthonormal matrix, or empty [] on first call.
    %           Existing orthonormal basis covering the first k columns of a.
    %
    %   r_in:   k-by-k upper triangular matrix, or empty [] on first call.
    %           Existing upper-triangular factor for the first k columns.
    %
    %   Output Parameters:
    %   ------------------
    %
    %   q_out:  m-by-n orthonormal matrix.
    %           Satisfies  a = q_out * r_out  (up to floating-point precision).
    %
    %   r_out:  n-by-n upper triangular matrix.
    %
    %   References:
    %   -----------
    %
    %   Adapted from qrupdate_gs.m in the GMRES-SDR reference implementation
    %   by the NLA Group, University of Manchester.
    %   https://github.com/nla-group/GMRES-SDR
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

    [n_rows, n_cols] = size(a);
    n_prev = size(q_in, 2);

    q_out = zeros(n_rows, n_cols);
    r_out = zeros(n_cols, n_cols);

    % q_in may have fewer rows than n_rows: when the Hessenberg gains one
    % row between consecutive calls (e.g. H grows from (j)x(j-1) to
    % (j+1)xj), the caller passes the grown matrix while q_in still has
    % the old row count.  Copy only the rows that exist; the new bottom
    % row of q_out stays zero, which is correct because the new Arnoldi
    % vector has zero projection onto the existing basis.
    n_q_rows = size(q_in, 1);
    q_out(1:n_q_rows, 1:n_prev) = q_in;
    r_out(1:n_prev, 1:n_prev) = r_in;

    for j = n_prev + 1:n_cols
        w = a(:, j);

        % Two passes of modified Gram-Schmidt for numerical stability.
        % The second pass corrects for cancellation in the first.
        for pass = 0:1 %#ok<*NODEF>
            for i = 1:j - 1
                proj = q_out(:, i)' * w;
                r_out(i, j) = r_out(i, j) + proj;
                w = w - q_out(:, i) * proj;
            end
        end

        r_out(j, j) = norm(w);
        q_out(:, j) = w / r_out(j, j);
    end

end
