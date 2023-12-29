function [HUpTri, g] = plane_rotations(H, beta)
    % Performs plane rotations.
    %
    %   Description:
    %   ------------
    %
    %   Implementation of plane rotations (a.k.a Givens rotation) following [1].
    %
    %   Syntaxis:
    %   ---------
    %
    %   [HUpTri, g] = plane_rotations(H, beta)
    %
    %   Input parameters:
    %   -----------------
    %
    %   H:          m+1-by-m matrix
    %               Upper Hessenberg matrix.
    %
    %   beta:       int
    %               Norm of the last cycle residual vector.
    %
    %   Output parameters:
    %   ------------------
    %
    %   HUpTri:     m+1-by-m matrix
    %               Upper-triangular Hessenberg matrix.
    %
    %   g:          m+1-by-1 vector
    %               Resulting right-hand-side vector.
    %
    %   Notes:
    %   ------
    %
    %   HUpTri and g (with the last rows deleted) will enter into the
    %   least-squares problem.
    %
    %   References:
    %   -----------
    %
    %   [1] Saad, Y. (2003). Iterative methods for sparse linear systems.
    %   Society for Industrial and Applied Mathematics.
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

    % Infer 'm' from the size of the upper Hessenber matrix H
    [~, m] = size(H);

    % Create rhs vector g
    g = zeros(m + 1, 1);
    g(1, 1) = beta;

    % Plane rotations
    for j = 1:m
        % Obtain sines and cosines
        s = H(j + 1, j) / (sqrt(abs(H(j + 1, j))^2 + H(j, j)^2));
        c = H(j, j) / (sqrt(abs(H(j + 1, j))^2 + H(j, j)^2));

        % Build rotation matrix
        P = eye(m + 1);
        P(j, j) = c;
        P(j + 1, j + 1) = c;
        P(j, j + 1) = s;
        P(j + 1, j) = -s;

        % Update HUpTri and g
        H = P * H;
        g = P * g;
    end

    HUpTri = H;

end
