function [H, v, W] = modified_gram_schmidt_arnoldi(A, v, m)
    % Modified Gram-Schmidt Arnoldi iteration
    %
    %   Description:
    %   ------------
    %
    %   TODO: Add a brief description...
    %
    %   Syntaxis:
    %   ---------
    %
    %   [H, W] = modified_gram_schmidt_arnoldi(A, v, m)
    %
    %   Input parameters:
    %   -----------------
    %
    %   A:  n-by-n matrix
    %       Matrix of coefficients.
    %
    %   v:  n-by-1 vector
    %       Normalized residual vector.
    %
    %   m:  int
    %       Restart parameter.
    %
    %   Output parameters:
    %   ------------------
    %
    %   H:  size?
    %       Upper Hessenberg matrix
    %
    %   W:  size?
    %       Orthonormal basis of Krylov subspace
    %
    %   v:  n-b1-1 vector
    %       Updated normalized residual vector.
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

    % Modified Gram Schmidt-Arnoldi
    for j = 1:m
        w(:, j) = A * v(:, j);
        for i = 1:j
            h(i, j) = w(:, j)' * v(:, i);
            w(:, j) = w(:, j) - h(i, j) * v(:, i);
        end
        h(j + 1, j) = norm(w(:, j));
        if h(j + 1, j) == 0
            m = j;
            h2 = zeros(m + 1, m); % Comment by JCC: VERIFY!!! (why?)
            for k = 1:m
                h2(:, k) = h(:, k);
            end
            h = h2;
        else
            v(:, j + 1) = w(:, j) / h(j + 1, j);
        end
    end

    H = h;
    W = w;

end
