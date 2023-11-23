function [H, V] = modified_gram_schmidt_arnoldi(A, v1, m)
    % Modified Gram-Schmidt Arnoldi iteration
    %
    %   Description:
    %   ------------
    %
    %   Implementation of the modified Gram-Schmidt Arnoldi iteration
    %   following [1], Algorithm 6.2.
    %
    %   Syntaxis:
    %   ---------
    %
    %   [H, V] = modified_gram_schmidt_arnoldi(A, v1, m)
    %
    %   Input parameters:
    %   -----------------
    %
    %   A:  n-by-n matrix
    %       Matrix of coefficients.
    %
    %   v1: n-by-1 vector
    %       Normalized residual vector.
    %
    %   m:  int
    %       Restart parameter.
    %
    %   Output parameters:
    %   ------------------
    %
    %   H:  m+1-by-m matrix
    %       Upper Hessenberg matrix
    %
    %   V:  n-by-m matrix
    %       Orthonormal basis of Krylov subspace
    %
    %   References:
    %   -----------
    %
    %   [1] Saad, Y. (2003). Iterative methods for sparse linear systems.
    %   Society for Industrial and Applied Mathematics.
    %
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

    % Initialize H and W
    [n, ~] = size(A);
    H = zeros(m+1, m);
    W = zeros(n, m);
    
    % Construct V matrix
    V = zeros(n, m + 1);
    V(:, 1) = v1;

    % Modified Gram Schmidt-Arnoldi
    for j = 1:m
        W(:, j) = A * V(:, j);
        for i = 1:j
            H(i, j) = W(:, j)' * V(:, i);
            W(:, j) = W(:, j) - H(i, j) * V(:, i);
        end
        H(j + 1, j) = norm(W(:, j));
        
        if H(j + 1, j) == 0
            % We have reached convergence. No need to continue.
            m = j;
            H = H(1:m + 1, 1:m);
            V = V(:, 1:m);
            return
        else
            V(:, j + 1) = W(:, j) / H(j + 1, j);
        end
    
    end

    % Slice matrix V since we are only interested in the first 'm' columns
    V = V(1:n, 1:m);
end
