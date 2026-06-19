% KrySBAS Octave benchmark - run by julia/benchmark.jl via subprocess.
% Each result line: RESULT,solver,matrix,n,nnz,time_s,iters,relres
%
% Usage (standalone):
%   octave --norc --no-gui matlab/benchmark.m
%
%   Copyright:
%   ----------
%
%   This file is part of the KrySBAS MATLAB Toolbox.
%
%   Copyright 2026 CC&MA - NIDTec - FP - UNA
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

warning('off', 'all');   % suppress near-singular triangular-solve warnings

script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, 'src')));
data_dir = fullfile(script_dir, '..', 'data');

matrices  = {'sherman1', 'sherman4', 'sherman5'};
m_krylov  = 100;
l_aug     = 3;
m_init    = 30;
m_stp     = 3;
tol       = 1e-9;
maxit     = 1000;
n_runs    = 3;

for mi = 1:length(matrices)
    mat_name = matrices{mi};
    mat_path = fullfile(data_dir, [mat_name '.mat']);
    if ~exist(mat_path, 'file')
        continue
    end
    load(mat_path, 'Problem');
    A = Problem.A;
    b = Problem.b;
    n_size = size(A, 1);
    nnz_A  = nnz(A);

    % --- lgmres ---
    times = zeros(n_runs, 1);
    for k = 1:n_runs
        [~, ~, rrv, ~, t] = lgmres(A, b, m_krylov, l_aug, tol, maxit);
        times(k) = t;
    end
    fprintf('RESULT,lgmres,%s,%d,%d,%.6f,%d,%.3e\n', ...
            mat_name, n_size, nnz_A, min(times), numel(rrv) - 1, rrv(end));

    % --- pd_gmres ---
    times = zeros(n_runs, 1);
    for k = 1:n_runs
        [~, ~, rrv, ~, t] = pd_gmres(A, b, m_init, [], m_stp, tol, maxit);
        times(k) = t;
    end
    fprintf('RESULT,pd_gmres,%s,%d,%d,%.6f,%d,%.3e\n', ...
            mat_name, n_size, nnz_A, min(times), numel(rrv) - 1, rrv(end));
end
