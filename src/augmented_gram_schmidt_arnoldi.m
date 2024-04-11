function [H, V, s] = ...
    augmented_gram_schmidt_arnoldi(A, v, m, appendV)
    % Modified Gram-Schmidt Arnoldi iteration, with 
    % Krylov subspace augmentation
    %
    %   Description:
    %   ------------
    %
    %   Implementation of the augmented Gram-Schmidt Arnoldi iteration
    %   following [1], Algorithm 2.1.
    %
    %   Syntaxis:
    %   ---------
    %
    %   [H, V, s] = ...
    %       augmented_gram_schmidt_arnoldi(A, v, m, appendV)
    %
    %   Input parameters:
    %   -----------------
    %
    %   A:          n-by-n matrix
    %               Matrix of coefficients.
    %
    %   v:          n-by-1 vector
    %               Normalized residual vector.
    %
    %   m:          int
    %               Restart parameter.
    %
    %   appendV:    n-by-k matrix
    %               Matrix of information vectors to append to the 
    %               current Krylov subspace
    %
    %   Output parameters:
    %   ------------------
    %
    %   H:          m+1-by-m matrix
    %               Upper Hessenberg matrix.
    %
    %   V:          n-by-m matrix
    %               Orthonormal basis of Krylov subspace.
    %
    %   s:          int
    %               Size of the new search subspace
    %
    %   References:
    %   -----------
    %
    %   [1] Saad, Y. (1997). Analysis of augmented Krylov subspace methods.
    %   SIAM Journal on Matrix Analysis and Applications, 18(2), 435-449.
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

    % Initialize matrices H and W
    [n, ~] = size(A);
    % For our augmented Krylov method, we compute the new 'size'
    % of the problem
    k = size(appendV, 2);
    s = m + k;
    
    H = zeros(s + 1, s);
    W = zeros(n, s);  %  W is only used internally by the algorithm

    % Construct the V matrix
    V = zeros(n, s + 1);
    V(:, 1) = v;  % this would be v1 if one follows Ref. [1].  
    
    % Modified Gram Schmidt-Arnoldi
    for j = 1:s
        if j <= m
            W(:, j) = A * V(:, j);
        else
            W(:, j) = A * appendV(:, k-(j-m-1));
        end

        for i = 1:j
            H(i, j) = W(:, j)' * V(:, i);
            W(:, j) = W(:, j) - H(i, j) * V(:, i);
        end
        H(j + 1, j) = norm(W(:, j));

        % From here on now, parameter should be 's'
        if H(j + 1, j) == 0
            % We have reached convergence. No need to continue.
            s = j;
            % Slice matrices and return outputs.
            H = H(1:s + 1, 1:s);
            V = V(:, 1:s);
            %mUpdated = s;
            return
        else
            V(:, j + 1) = W(:, j) / H(j + 1, j);
        end

    end

    % Slice matrix V since we are only interested in the first 's' columns
    V = V(:, 1:s);

end
