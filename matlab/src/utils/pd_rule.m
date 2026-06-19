function miter = pd_rule(m, n, mInitial, mMin, mMax, mStep,  ...
                         res, iter, alphaP, alphaD)
    % Rule for the Proportional-Derivative law. Algorithm 1 of [1].
    %
    %   Signature:
    %   ----------
    %
    %   miter = (m, n, mInitial, mMin, mMax, mStep,  ...
    %                     res, iter, alphaP, alphaD)
    %
    %   Input Parameters:
    %   -----------------
    %
    %   m:          int
    %               Restart parameter at the last cycle
    %
    %   n:          int
    %               size of the square matrix A
    %
    %   mInitial:   int
    %               Initial restart parameter at the last cycle.
    %
    %   mMin:       int
    %               Minimum value of the restart parameter m.
    %
    %   res         number of cycles-by-1 vector
    %               Vector of residual norms of every outer iteration
    %               (cycles).
    %
    %   iter        number of cycles-by-1 vector
    %               Number of restart cycles previous.
    %
    %   mStep:      int
    %               Step size for increasing the mInitial when m < mMinMax(1).
    %
    %   mMax:       int
    %               Maximum value of the restart paramter m.
    %
    %   maxit:      int
    %               Maximum number of outer iterations.
    %
    %   alphaP:     int
    %               Proportional coefficient from PD controller.
    %
    %   alphaD:     int
    %               Derivative coefficient from PD controller.
    %
    %
    %   Output parameters:
    %   ------------------
    %
    %   miter:      1-by-2 vector
    %               Restart parameter and Initial restart parameter at the
    %               new cycle.
    %
    %    %   References:
    %   -----------
    %
    %   [1] Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A
    %   proportional-derivative control strategy for restarting the GMRES(m)
    %   algorithm. Journal of Computational and Applied Mathematics,
    %   337, 209-224.
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

    if iter > 3

        mj = m + ceil( ...
                      alphaP * (res(iter) / res(iter - 1)) + ...
                      alphaD * ( ...
                                (res(iter) - res(iter - 2)) / ...
                                (2 * res(iter - 1)) ...
                               ) ...
                     );

    elseif iter > 2
        mj = m + ceil(alphaP * (res(iter) / res(iter - 1)));

    else
        mj = mInitial;
    end

    if mj < mMin
        mInitial = mInitial + mStep;
        mj = mInitial;
    end

    if mj > mMax
        mj = mMax;
    end

    if mInitial + mStep > n
        mj = n;
    end

    miter = [mj mInitial];
