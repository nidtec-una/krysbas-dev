% compare_solvers.m
%
%   Description:
%   ------------
%
%   Compares the four KrySBAS solvers
%
%       GMRES-E(m, d)  -- restarted GMRES augmented with d harmonic Ritz
%                         vectors (Morgan, 1995)
%       LGMRES(m, l)   -- restarted GMRES augmented with l error
%                         approximation vectors (Baker et al., 2005)
%       PD-GMRES(m)    -- restarted GMRES with a PD controller that adapts
%                         m automatically (Nunez et al., 2018)
%       GMRES-DR(m, k) -- restarted GMRES with deflated restarting using
%                         k harmonic Ritz vectors (Morgan, 2002)
%
%   For each test matrix the script produces one figure showing the base-10
%   logarithm of the relative residual norm versus the restart cycle index
%   for all four solvers.  A horizontal dashed line marks the convergence
%   tolerance.
%
%   Usage:
%   ------
%
%   Run from the matlab/ directory (or from anywhere, the script locates
%   its own paths automatically):
%
%       >> compare_solvers
%
%   To change matrices or solver parameters, edit the CONFIGURATION section
%   below.
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
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%
%   KrySBAS is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
%   for more details.
%
%   You should have received a copy of the GNU General Public License along
%   with this file.  If not, see <http://www.gnu.org/licenses/>.

% =========================================================================
% ----> Self-contained path setup
%
% Locate the solver source and shared data directory relative to this
% script so that compare_solvers.m works regardless of the MATLAB current
% directory.
% =========================================================================

script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, 'src')));
data_dir = fullfile(script_dir, '..', 'data');

% =========================================================================
% ----> CONFIGURATION
%
% Edit this section to change which matrices are loaded and which solver
% parameters are used.  All solvers use the same m, tolerance, and maximum
% number of cycles so that the comparison is fair.
% =========================================================================

% Test matrices to compare (must exist as <name>.mat inside data/).
% Each .mat file must contain a struct 'Problem' with fields A and b.
matrices = {'embree3', 'sherman1', 'sherman4'};

% Shared restart dimension and convergence tolerance
m    = 27;   % Krylov subspace dimension (restart parameter)
tol  = 1e-9; % relative residual tolerance
maxit = 500; % maximum number of restart cycles

% GMRES-E: number of harmonic Ritz vectors to append (augmented dimension)
d_gmrese = 3;

% LGMRES: number of error approximation vectors to append
l_lgmres = 3;

% PD-GMRES: initial restart parameter and step size for the PD controller
m_init_pd = m;
m_step_pd = 3;

% GMRES-DR: number of harmonic Ritz vectors to deflate and recycle
k_gmresdr = 3;

% =========================================================================
% ----> Colour and line-style scheme (one style per solver)
% =========================================================================

styles = struct( ...
                'gmres_e', struct('color', [0.00 0.45 0.70], ...
                                  'line', '-', 'marker', 'o', ...
                                  'label', 'GMRES-E'), ...
                'lgmres', struct('color', [0.84 0.37 0.00], ...
                                 'line', '--', 'marker', 's', ...
                                 'label', 'LGMRES'), ...
                'pd_gmres', struct('color', [0.00 0.62 0.45], ...
                                   'line', ':', 'marker', '^', ...
                                   'label', 'PD-GMRES'), ...
                'gmres_dr', struct('color', [0.80 0.47 0.65], ...
                                   'line', '-.', 'marker', 'd', ...
                                   'label', 'GMRES-DR') ...
               );

% =========================================================================
% ----> Main loop: one figure per matrix
% =========================================================================

for mi = 1:length(matrices)

    mat_name = matrices{mi};
    mat_path = fullfile(data_dir, [mat_name '.mat']);

    % Check that the file exists before trying to load it
    if ~exist(mat_path, 'file')
        warning('Matrix file not found, skipping: %s', mat_path);
        continue
    end

    load(mat_path, 'Problem');
    A = Problem.A;
    b = Problem.b;
    n = size(A, 1);

    fprintf('\n--- %s  (n = %d, nnz = %d) ---\n', mat_name, n, nnz(A));

    % ------------------------------------------------------------------
    % Run each solver and collect (relresvec, flag, time).
    % ------------------------------------------------------------------

    % GMRES-E(m, d)
    [~, flag_e, rrv_e, ~, t_e] = gmres_e(A, b, m, d_gmrese, tol, maxit);
    fprintf('  GMRES-E   : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_e, numel(rrv_e) - 1, t_e);

    % LGMRES(m, l)
    [~, flag_l, rrv_l, ~, t_l] = lgmres(A, b, m, l_lgmres, tol, maxit);
    fprintf('  LGMRES    : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_l, numel(rrv_l) - 1, t_l);

    % PD-GMRES
    [~, flag_p, rrv_p, ~, t_p] = pd_gmres(A, b, m_init_pd, [], ...
                                          m_step_pd, tol, maxit);
    fprintf('  PD-GMRES  : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_p, numel(rrv_p) - 1, t_p);

    % GMRES-DR(m, k)
    [~, flag_dr, rrv_dr, ~, t_dr] = gmres_dr(A, b, m, k_gmresdr, ...
                                             tol, maxit);
    fprintf('  GMRES-DR  : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_dr, numel(rrv_dr) - 1, t_dr);

    % ------------------------------------------------------------------
    % Build cycle-index vectors (0-based: entry 1 = before iteration 1)
    % ------------------------------------------------------------------
    cycles_e  = 0:numel(rrv_e)  - 1;
    cycles_l  = 0:numel(rrv_l)  - 1;
    cycles_p  = 0:numel(rrv_p)  - 1;
    cycles_dr = 0:numel(rrv_dr) - 1;

    % ------------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------------
    figure('Name', ['Solver comparison -- ' mat_name], ...
           'NumberTitle', 'off');
    hold on;

    semilogy(cycles_e,  rrv_e, ...
             'Color',     styles.gmres_e.color,  ...
             'LineStyle', styles.gmres_e.line,   ...
             'Marker',    styles.gmres_e.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,  ...
             'DisplayName', styles.gmres_e.label);

    semilogy(cycles_l,  rrv_l, ...
             'Color',     styles.lgmres.color,   ...
             'LineStyle', styles.lgmres.line,    ...
             'Marker',    styles.lgmres.marker,  ...
             'LineWidth', 1.5, 'MarkerSize', 5,  ...
             'DisplayName', styles.lgmres.label);

    semilogy(cycles_p,  rrv_p, ...
             'Color',     styles.pd_gmres.color,  ...
             'LineStyle', styles.pd_gmres.line,   ...
             'Marker',    styles.pd_gmres.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,   ...
             'DisplayName', styles.pd_gmres.label);

    semilogy(cycles_dr, rrv_dr, ...
             'Color',     styles.gmres_dr.color,  ...
             'LineStyle', styles.gmres_dr.line,   ...
             'Marker',    styles.gmres_dr.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,   ...
             'DisplayName', styles.gmres_dr.label);

    % Horizontal tolerance line
    x_max = max([cycles_e(end), cycles_l(end), ...
                 cycles_p(end), cycles_dr(end)]);
    semilogy([0, x_max], [tol, tol], ...
             'k--', 'LineWidth', 1.0, ...
             'DisplayName', sprintf('tol = %.0e', tol));

    hold off;

    xlabel('Restart cycle', 'FontSize', 12);
    ylabel('Relative residual norm', 'FontSize', 12);
    title(sprintf('Solver comparison  --  %s  (n = %d)', mat_name, n), ...
          'FontSize', 13);
    legend('Location', 'northeast', 'FontSize', 10);
    grid on;
    set(gca, 'YScale', 'log');

end
