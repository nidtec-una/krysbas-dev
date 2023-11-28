function miter = pd_rule(m, minitial, mmin, ...
                         res, iter, mstep, mmax, alpha, delta)
    % Rule for the Proportional-Derivative law.
    %
    %   TODO: Complete rest of the docstring.
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
                      alpha * (res(iter, :) / res(iter - 1, :)) + ...
                      delta * ( ...
                               (res(iter, :) - res(iter - 2, :)) / ...
                               (2 * res(iter - 1, :)) ...
                              ) ...
                     );


    elseif iter > 2
        mj = m + ceil(alpha * (res(iter, :) / res(iter - 1, :)));
        

    else
        mj = minitial;
    end

    % ap=res(iter,:)/res(iter-1,:)
    % if iter >3
    % ad=(res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))
    % end

    % if mj <= 0
    %      mj=abs(mj) + mstep;
    % end
    if mj < mmin
        minitial = minitial + mstep;
        mj = minitial;
    end

    if mj > mmax
        %    miter=minitial;
        %     miter=miter-mstep;
        % mj=mmax-mstep;
        mj = mmax;
    end

    miter = [mj minitial];
    
