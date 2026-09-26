%% Problem 1: high-frequency benchmark of a solid copper wire (MATLAB side)
% Runs standalone. If the COMSOL exports (radial_200kHz.csv,
% mesh_convergence.csv, frequency_sweep.csv, results_summary.json) sit in
% this same folder, the script overlays them and prints the comparison.
% LiveLink is not required.

clear; clc; close all;
thisDir = fileparts(mfilename('fullfile'));
if ~isempty(thisDir)
    cd(thisDir);
end

%% Common parameters
sigma = 5.8e7;                 % Copper conductivity, S/m
rho = 1 / sigma;               % Resistivity, ohm*m
mu0 = 4 * pi * 1e-7;           % Vacuum permeability, H/m
mur = 1;
mu = mu0 * mur;
f = 200e3;                     % Frequency, Hz
omega = 2 * pi * f;
I_rms = 20;                    % RMS current, A
D = 2e-3;                      % Wire diameter, m
a = D / 2;                     % Wire radius, m
L = 1;                         % Evaluation length, m

%% Analytical field solution (exp(+j*omega*t), RMS phasors)
% J(r) = C*J0(k*r), k = (1-j)/delta, C set by total RMS current.
delta = sqrt(2 / (omega * mu * sigma));
k = (1 - 1i) / delta;
Cconst = I_rms * k / (2 * pi * a * besselj(1, k * a));
nr = 4001;
r = linspace(0, a, nr).';
J = Cconst * besselj(0, k * r);
Jmag = abs(J);

%% Task 1: direct results via disp
Rdc = L / (sigma * pi * a^2);
Jdc = I_rms / (pi * a^2);
Zint = L * k / (2 * pi * a * sigma) * besselj(0, k * a) / besselj(1, k * a);
Rac_exact = real(Zint);

disp('=== Task 1: direct results ===');
disp(['Frequency f = ', num2str(f / 1e3, '%.6g'), ' kHz']);
disp(['Skin depth delta = ', num2str(delta * 1e3, '%.9g'), ' mm']);
disp(['Wire radius a = ', num2str(a * 1e3, '%.6g'), ' mm']);
disp(['a/delta = ', num2str(a / delta, '%.9g')]);
disp(['DC resistance Rdc = ', num2str(Rdc * 1e3, '%.9g'), ' mOhm/m']);
disp(['DC current density Jdc = ', num2str(Jdc / 1e6, '%.9g'), ' A/mm^2']);
disp(['Exact AC resistance (Bessel surface impedance) = ', ...
    num2str(Rac_exact * 1e3, '%.9g'), ' mOhm/m']);

%% Task 2(a): radial current-density plot (one figure)
haveComsolProfile = isfile('radial_200kHz.csv');
fig = figure('Color', 'w', 'Name', 'Problem 1 - radial current density');
plot(r * 1e3, Jmag / 1e6, 'LineWidth', 1.6); hold on;
if haveComsolProfile
    comsolData = readmatrix('radial_200kHz.csv');
    plot(comsolData(:, 1) * 1e3, comsolData(:, 3) / 1e6, '--', 'LineWidth', 1.2);
    legend('MATLAB |J(r)| (analytic)', 'COMSOL FEM |J(r)| (radial_200kHz.csv)', ...
        'Location', 'northeast');
else
    legend('MATLAB |J(r)| (analytic)', 'Location', 'northeast');
end
grid on; box on;
xlabel('Radius r (mm)');
ylabel('|J(r)| (A/mm^2)');
title('Solid copper wire: RMS current-density magnitude, 200 kHz');
saveas(fig, fullfile(pwd, 'q1_matlab_vs_comsol_radial.png'));

%% Task 2(b): COMSOL data interface (no LiveLink: CSV exchange)
writematrix([r, real(J), imag(J), Jmag, angle(J)], ...
    fullfile(pwd, 'matlab_theory_radial_profile.csv'));
writematrix([f, sigma, mu, I_rms, D, L, delta, Rdc, Rac_exact], ...
    fullfile(pwd, 'matlab_parameters.csv'));

disp('=== Task 2: plot and COMSOL data interface ===');
disp('Saved q1_matlab_vs_comsol_radial.png,');
disp('matlab_theory_radial_profile.csv and matlab_parameters.csv.');
if haveComsolProfile
    JrmsFEM = interp1(comsolData(:, 1), comsolData(:, 3), r, 'linear', 'extrap');
    JzFEM = interp1(comsolData(:, 1), ...
        comsolData(:, 4) + 1i * comsolData(:, 5), r, 'linear', 'extrap');
    errMag = max(abs(JrmsFEM - Jmag)) / Jmag(end) * 100;
    errCplx = max(abs(JzFEM - J)) / Jmag(end) * 100;
    fprintf('COMSOL vs MATLAB |J| max error: %.6f%% of surface value\n', errMag);
    fprintf('COMSOL vs MATLAB complex Jz max error: %.6f%% of surface value\n', errCplx);
else
    disp('radial_200kHz.csv not found: comparison skipped.');
end

%% Task 3(a): MATLAB loss calculation only
% RMS phasors: P = (1/sigma)*integral(|J|^2 dV), dV = 2*pi*r*dr*L.
P_ac = (L / sigma) * 2 * pi * trapz(r, Jmag.^2 .* r);
Rac = P_ac / I_rms^2;
ratio = Rac / Rdc;

disp('=== Task 3: MATLAB loss calculation ===');
disp(['AC loss P_ac = ', num2str(P_ac, '%.9g'), ' W/m']);
disp(['AC resistance Rac (trapz, Nr=4001) = ', num2str(Rac * 1e3, '%.9g'), ' mOhm/m']);
disp(['Rac/Rdc = ', num2str(ratio, '%.9g')]);

%% Task 3(b): summary export; contour plots are produced by COMSOL
writematrix([P_ac, Rac, Rdc, ratio, Rac_exact], ...
    fullfile(pwd, 'matlab_summary.csv'));
disp('Saved matlab_summary.csv; cross-section contour comes from COMSOL.');
if isfile('comsol_cross_section.png')
    disp('COMSOL contour available: comsol_cross_section.png');
end

%% Task 4(a): LaTeX formulas printed directly
disp('=== Task 4: LaTeX formulas and numerical solution ===');
disp('delta = sqrt(2/(omega*mu*sigma)) = sqrt(rho/(pi*f*mu))');
disp('J(r) = C*besselj(0,k*r),  k = (1-j)/delta');
disp('C = I*k/(2*pi*a*besselj(1,k*a))');
disp('P_ac = (L/sigma)*2*pi*integral_0^a |J(r)|^2*r dr');
disp('R_ac = P_ac/I^2,  R_dc = L/(sigma*pi*a^2)');
disp('Z_int = L*k/(2*pi*a*sigma) * besselj(0,k*a)/besselj(1,k*a),  R_ac = Re(Z_int)');

%% Task 4(b): numerical solution plus convergence and COMSOL comparison
grid_counts = [101, 201, 401, 801, 1601, 3201, 6401];
convergence = zeros(numel(grid_counts), 4);
for q = 1:numel(grid_counts)
    rq = linspace(0, a, grid_counts(q)).';
    Jq = Cconst * besselj(0, k * rq);
    Pq = (L / sigma) * 2 * pi * trapz(rq, abs(Jq).^2 .* rq);
    Raq = Pq / I_rms^2;
    convergence(q, :) = [grid_counts(q), Pq, Raq, Raq / Rdc];
end
writematrix(convergence, fullfile(pwd, 'matlab_grid_convergence.csv'));
disp('MATLAB convergence (Nr, P_ac, R_ac, R_ac/R_dc):');
disp(convergence);

if isfile('mesh_convergence.csv')
    meshData = readmatrix('mesh_convergence.csv');
    fprintf('COMSOL final mesh: Rac = %.9f mOhm/m, Rac/Rdc = %.9f\n', ...
        meshData(end, 8) * 1e3, meshData(end, 9));
    fprintf('Rac difference MATLAB(trapz) - COMSOL(FEM): %.6e mOhm/m\n', ...
        (Rac - meshData(end, 8)) * 1e3);
end
if isfile('results_summary.json')
    summaryJSON = jsondecode(fileread('results_summary.json'));
    fprintf('COMSOL summary: delta = %.9g um, Rac_FEM = %.9g mOhm, ratio = %.9g\n', ...
        summaryJSON.delta_um, summaryJSON.Rac_FEM_mohm, summaryJSON.ratio_FEM);
    if isfield(summaryJSON, 'skin_layer_efold_FEM_um')
        fprintf('1/e skin layer: FEM = %.4f um, Bessel = %.4f um, plane delta = %.4f um\n', ...
            summaryJSON.skin_layer_efold_FEM_um, ...
            summaryJSON.skin_layer_efold_Bessel_um, summaryJSON.delta_um);
    end
end

disp('=== Final numerical solution ===');
disp(['delta = ', num2str(delta, '%.12g'), ' m']);
disp(['P_ac = ', num2str(P_ac, '%.12g'), ' W/m']);
disp(['R_ac = ', num2str(Rac, '%.12g'), ' ohm/m (Bessel exact: ', ...
    num2str(Rac_exact, '%.12g'), ')']);
disp(['R_ac/R_dc = ', num2str(ratio, '%.12g')]);
