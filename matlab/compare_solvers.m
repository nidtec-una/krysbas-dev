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
%   Copyright 2026 CC&MA - NIDTec - FP - UNA
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
matrices = {'sherman4'};

% Shared restart dimension and convergence tolerance
m    = 27;   % Krylov subspace dimension (restart parameter)
tol  = 1e-9; % relative residual tolerance
maxit = 500; % maximum number of restart cycles

% GMRES-E: number of harmonic Ritz vectors to append (augmented dimension)
d_gmrese = 3;

% LGMRES: number of error approximation vectors to append
l_lgmres = 3;

% PD-GMRES: initial restart parameter, step size, and PD coefficients.
% alpha_p (proportional) and alpha_d (derivative) default to [-3, 5]
% following Nunez et al. (2018), p. 217.
m_init_pd = m;
m_step_pd = 3;
alpha_p   = -3;
alpha_d   =  5;

% GMRES-DR: number of harmonic Ritz vectors to deflate and recycle
k_gmresdr = 3;

% =========================================================================
% ----> Colour and line-style scheme (one style per solver)
% =========================================================================

styles = struct( ...
                'gmres_m', struct('color', [0.60 0.60 0.60], ...
                                  'line', '-', 'marker', 'x', ...
                                  'label', 'GMRES(m)'), ...
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

    % Plain restarted GMRES(m) -- built-in, used as baseline.
    % gmres() returns resvec with one entry per inner iteration; we
    % downsample to one entry per restart cycle by stepping every m
    % positions.  If convergence occurs mid-cycle the final entry is
    % appended so the last point is always included.
    tic_gm = tic();
    [~, flag_gm, ~, ~, resvec_gm] = gmres(A, b, m, tol, maxit);
    t_gm = toc(tic_gm);
    idx_gm = unique([1:m:numel(resvec_gm), numel(resvec_gm)]);
    rrv_gm = resvec_gm(idx_gm);
    fprintf('  GMRES(m)  : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_gm, numel(rrv_gm) - 1, t_gm);

    % GMRES-E(m, d)
    [~, flag_e, rrv_e, ~, t_e] = gmres_e(A, b, m, d_gmrese, tol, maxit);
    fprintf('  GMRES-E   : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_e, numel(rrv_e) - 1, t_e);

    % LGMRES(m, l)
    [~, flag_l, rrv_l, ~, t_l] = lgmres(A, b, m, l_lgmres, tol, maxit);
    fprintf('  LGMRES    : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_l, numel(rrv_l) - 1, t_l);

    % PD-GMRES (alphaPD passed explicitly so the legend shows actual values)
    [~, flag_p, rrv_p, ~, t_p] = pd_gmres(A, b, m_init_pd, [], ...
                                          m_step_pd, tol, maxit, ...
                                          [], [alpha_p; alpha_d]);
    fprintf('  PD-GMRES  : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_p, numel(rrv_p) - 1, t_p);

    % GMRES-DR(m, k)
    [~, flag_dr, rrv_dr, ~, t_dr] = gmres_dr(A, b, m, k_gmresdr, ...
                                             tol, maxit);
    fprintf('  GMRES-DR  : flag=%d  cycles=%3d  time=%.3fs\n', ...
            flag_dr, numel(rrv_dr) - 1, t_dr);

    % ------------------------------------------------------------------
    % Build legend labels with solver parameters and CPU time.
    % LaTeX interpreter is used so \alpha_P, \alpha_D render as Greek.
    % ------------------------------------------------------------------
    lbl_gm = ['$\mathrm{GMRES}(m)\ ' ...
              '(m=' num2str(m) ...
              ',\ t=' sprintf('%.2f', t_gm) '\ \mathrm{s})$'];

    lbl_e  = ['$\mathrm{GMRES\mbox{-}E}\ ' ...
              '(m=' num2str(m) ',\ d=' num2str(d_gmrese) ...
              ',\ t=' sprintf('%.2f', t_e) '\ \mathrm{s})$'];

    lbl_l  = ['$\mathrm{LGMRES}\ ' ...
              '(m=' num2str(m) ',\ l=' num2str(l_lgmres) ...
              ',\ t=' sprintf('%.2f', t_l) '\ \mathrm{s})$'];

    lbl_p  = ['$\mathrm{PD\mbox{-}GMRES}\ ' ...
              '(m=' num2str(m_init_pd) ...
              ',\ \alpha_P=' num2str(alpha_p) ...
              ',\ \alpha_D=' num2str(alpha_d) ...
              ',\ t=' sprintf('%.2f', t_p) '\ \mathrm{s})$'];

    lbl_dr = ['$\mathrm{GMRES\mbox{-}DR}\ ' ...
              '(m=' num2str(m) ',\ k=' num2str(k_gmresdr) ...
              ',\ t=' sprintf('%.2f', t_dr) '\ \mathrm{s})$'];

    lbl_tol = ['$\mathrm{tol} = ' sprintf('%g', tol) '$'];

    % ------------------------------------------------------------------
    % Build cycle-index vectors (0-based: entry 1 = before iteration 1)
    % ------------------------------------------------------------------
    cycles_gm = 0:numel(rrv_gm) - 1;
    cycles_e  = 0:numel(rrv_e)  - 1;
    cycles_l  = 0:numel(rrv_l)  - 1;
    cycles_p  = 0:numel(rrv_p)  - 1;
    cycles_dr = 0:numel(rrv_dr) - 1;

    % ------------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------------
    figure('Name', ['Solver comparison -- ' mat_name], ...
           'NumberTitle', 'off', ...
           'Units', 'pixels', 'Position', [50, 50, 1150, 520]);
    hold on;

    semilogy(cycles_gm, rrv_gm, ...
             'Color',     styles.gmres_m.color,  ...
             'LineStyle', styles.gmres_m.line,   ...
             'Marker',    styles.gmres_m.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,  ...
             'DisplayName', lbl_gm);

    semilogy(cycles_e,  rrv_e, ...
             'Color',     styles.gmres_e.color,  ...
             'LineStyle', styles.gmres_e.line,   ...
             'Marker',    styles.gmres_e.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,  ...
             'DisplayName', lbl_e);

    semilogy(cycles_l,  rrv_l, ...
             'Color',     styles.lgmres.color,   ...
             'LineStyle', styles.lgmres.line,    ...
             'Marker',    styles.lgmres.marker,  ...
             'LineWidth', 1.5, 'MarkerSize', 5,  ...
             'DisplayName', lbl_l);

    semilogy(cycles_p,  rrv_p, ...
             'Color',     styles.pd_gmres.color,  ...
             'LineStyle', styles.pd_gmres.line,   ...
             'Marker',    styles.pd_gmres.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,   ...
             'DisplayName', lbl_p);

    semilogy(cycles_dr, rrv_dr, ...
             'Color',     styles.gmres_dr.color,  ...
             'LineStyle', styles.gmres_dr.line,   ...
             'Marker',    styles.gmres_dr.marker, ...
             'LineWidth', 1.5, 'MarkerSize', 5,   ...
             'DisplayName', lbl_dr);

    % Horizontal tolerance line
    x_max = max([cycles_gm(end), cycles_e(end), ...
                 cycles_l(end), cycles_p(end), cycles_dr(end)]);
    semilogy([0, x_max], [tol, tol], ...
             'k--', 'LineWidth', 1.0, 'DisplayName', lbl_tol);

    hold off;

    xlabel('Restart cycle', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('$\|r\| / \|r_0\|$', 'FontSize', 12, 'Interpreter', 'latex');
    title(['Solver comparison -- ' mat_name ...
           ' $(n = ' num2str(n) ')$'], ...
          'FontSize', 13, 'Interpreter', 'latex');
    legend('Location', 'eastoutside', 'FontSize', 10, ...
           'Interpreter', 'latex');
    grid on;
    set(gca, 'YScale', 'log');

end
