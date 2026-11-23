% BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
% Teacher in charge: Marko Matikainen (MKM)
%
% Mesh convergence of eigenfrequencies: relative error between successive meshes
%
%   e_i = |omega_i - omega_{i-1}| / omega_i * 100 %
%
%   omega_i     : eigenfrequency of mesh i
%   omega_{i-1} : eigenfrequency of the previous (coarser) mesh
%
% Example: free plane structure (L = 2 m, H = 0.4 m, Lz = 0.4 m, E = 1 GPa,
% nu = 0.3, rho = 7850 kg/m^3), Q4 elements, consistent mass matrix.
% Meshes 4x1, 8x2, 16x4, 32x8 (ScaleMesh = 1, 2, 4, 8).
% Full Gauss integration (KlocQ4xiFullGauss) vs. SRI (KlocQ4psR).
%
% The eigenfrequencies can be
%   1) typed in below (ReadFromFile = false), or
%   2) read from a csv file (ReadFromFile = true), see omega_results_Q4.csv
%
% NOTE: compare modes of the same type! The mode order may change between
% meshes (e.g. 4x1: axial mode before bending mode 2). Check the mode shapes.

clear; close all; clc;

ReadFromFile = false;

NumElems  = [4 16 64 256];                    % number of elements
MeshNames = {'4x1','8x2','16x4','32x8'};
ModeNames = {'Bending 1','Bending 2','Bending 3','Axial 1','Axial 2'};

if ~ReadFromFile
    % Eigenfrequencies omega [rad/s]; rows: modes, columns: meshes 4x1 ... 32x8
    omegaFull = [ 281.93   226.14   210.03   205.80;     % Bending 1
                  736.78   554.82   499.02   484.31;     % Bending 2
                 1273.40   974.98   845.68   811.31;     % Bending 3
                  575.43   563.77   560.79   560.03;     % Axial 1
                 1236.99  1145.84  1121.65  1115.53];    % Axial 2

    omegaSRI  = [ 231.42   212.35   206.45   204.90;     % Bending 1
                  634.95   526.42   491.48   482.38;     % Bending 2
                 1189.28   934.56   834.13   808.30;     % Bending 3
                  575.26   563.70   560.77   560.03;     % Axial 1
                 1235.73  1145.13  1121.43  1115.47];    % Axial 2
else
    % csv: first row = header, columns: Mode, Integration, 4x1, 8x2, 16x4, 32x8
    T = readtable('omega_results_Q4.csv');
    omegaFull = T{strcmp(T.Integration,'full'), 3:end};
    omegaSRI  = T{strcmp(T.Integration,'SRI'),  3:end};
end

ns    = numel(NumElems);
nmode = numel(ModeNames);

%% Relative error between successive meshes [%]
eFull = abs(omegaFull(:,2:ns) - omegaFull(:,1:ns-1)) ./ omegaFull(:,2:ns) * 100;
eSRI  = abs(omegaSRI(:,2:ns)  - omegaSRI(:,1:ns-1))  ./ omegaSRI(:,2:ns)  * 100;

%% Table
fprintf('                          omega [rad/s]                         |  e_i [%%]\n');
fprintf('Mode       Int.        4x1       8x2      16x4      32x8 |    8x2   16x4   32x8\n');
for k = 1:nmode
    fprintf('%-10s full %9.2f %9.2f %9.2f %9.2f | %6.2f %6.2f %6.2f\n', ...
        ModeNames{k}, omegaFull(k,:), eFull(k,:));
    fprintf('%-10s SRI  %9.2f %9.2f %9.2f %9.2f | %6.2f %6.2f %6.2f\n', ...
        '', omegaSRI(k,:), eSRI(k,:));
end

% The same as a MATLAB table
Results = table(repelem(ModeNames',2), repmat({'full';'SRI'},nmode,1), ...
    reshape([omegaFull(:,ns) omegaSRI(:,ns)]',[],1), ...
    reshape([eFull(:,end) eSRI(:,end)]',[],1), ...
    'VariableNames', {'Mode','Integration','omega_32x8','e_32x8_percent'})

%% Convergence plot (log-log): full Gauss solid, SRI dashed
x       = NumElems(2:ns);                     % finer mesh of each pair
colors  = lines(nmode);
markers = {'o','s','^','d','v'};

figure; hold on
for k = 1:nmode
    plot(x, eFull(k,:), '-',  'Color', colors(k,:), 'Marker', markers{k}, ...
        'MarkerFaceColor', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', [ModeNames{k} ', full Gauss']);
    plot(x, eSRI(k,:),  '--', 'Color', colors(k,:), 'Marker', markers{k}, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
        'DisplayName', [ModeNames{k} ', SRI']);
end
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on; box on
xticks(x); xticklabels({'16 (8x2)','64 (16x4)','256 (32x8)'});
xlabel('Number of elements (finer mesh)');
ylabel('Relative error e_i [%]');
title('e_i = |\omega_i - \omega_{i-1}| / \omega_i \cdot 100 %');
legend('Location', 'eastoutside');
