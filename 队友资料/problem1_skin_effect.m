%% Problem 1: high-frequency simulation of a solid copper wire
% The script follows the four tasks in the problem statement.
% It uses the analytical cylindrical-conductor solution as the MATLAB
% reference model and exports the radial profile for COMSOL comparison.

clear; clc; close all;

%% Common parameters
sigma = 5.8e7;                 % Copper conductivity, S/m
rho = 1 / sigma;               % Resistivity, ohm*m
mu0 = 4 * pi * 1e-7;          % Vacuum permeability, H/m
mur = 1;
mu = mu0 * mur;
f = 200e3;                     % Frequency, Hz
omega = 2 * pi * f;
I_rms = 20;                    % RMS current, A
D = 2e-3;                      % Wire diameter, m
a = D / 2;                     % Wire radius, m
L = 1;                         % Evaluation length, m

%% Task 1: skin depth and direct result output
delta = sqrt(2 / (omega * mu * sigma));
Rdc = L / (sigma * pi * a^2);
Jdc = I_rms / (pi * a^2);

disp('=== Task 1: direct results ===');
disp(['Frequency f = ', num2str(f / 1e3, '%.6g'), ' kHz']);
disp(['Skin depth delta = ', num2str(delta * 1e3, '%.9g'), ' mm']);
disp(['Wire radius a = ', num2str(a * 1e3, '%.6g'), ' mm']);
disp(['a/delta = ', num2str(a / delta, '%.9g')]);
disp(['DC resistance Rdc = ', num2str(Rdc * 1e3, '%.9g'), ' mOhm']);
disp(['DC current density Jdc = ', num2str(Jdc / 1e6, '%.9g'), ' A/mm^2']);

%% Task 2(a): radial current-density plot
% With exp(+j*omega*t), the diffusion wavenumber is k=(1-j)/delta.
% Reversing the time convention conjugates J but leaves |J| and loss unchanged.
k = (1 - 1i) / delta;
nr = 4001;
r = linspace(0, a, nr).';

% J(r) = C*J0(k*r), normalized to the specified RMS total current.
C = I_rms * k / (2 * pi * a * besselj(1, k * a));
J = C * besselj(0, k * r);
Jmag = abs(J);

figure('Color', 'w', 'Name', 'Problem 1 - radial current density');
plot(r * 1e3, Jmag / 1e6, 'LineWidth', 1.6);
grid on; box on;
xlabel('Radius r (mm)');
ylabel('|J(r)| (A/mm^2)');
title('Solid copper wire: RMS current-density magnitude');

%% Task 2(b): COMSOL LiveLink-compatible exports
% These CSV files can also be read by a COMSOL model or by LiveLink for MATLAB.
profile = [r, real(J), imag(J), Jmag, angle(J)];
parameters = [f, sigma, mu, I_rms, D, L, delta, Rdc];
writematrix(profile, fullfile(pwd, 'problem1_radial_profile.csv'));
writematrix(parameters, fullfile(pwd, 'problem1_parameters.csv'));

disp('=== Task 2: plot and COMSOL interface ===');
disp('Saved problem1_radial_profile.csv and problem1_parameters.csv.');
disp('The files can be imported into COMSOL, or read through LiveLink if licensed.');

%% Task 3(a): MATLAB loss calculation only
% RMS phasors are used, so P = (1/sigma)*integral(|J|^2 dV).
% For a unit-length wire, dV = 2*pi*r*dr*L.
P_ac = (L / sigma) * 2 * pi * trapz(r, Jmag.^2 .* r);
Rac = P_ac / I_rms^2;
ratio = Rac / Rdc;

disp('=== Task 3: MATLAB loss calculation ===');
disp(['AC loss P_ac = ', num2str(P_ac, '%.9g'), ' W/m']);
disp(['AC resistance Rac = ', num2str(Rac * 1e3, '%.9g'), ' mOhm/m']);
disp(['Rac/Rdc = ', num2str(ratio, '%.9g')]);

%% Task 3(b): export values for COMSOL visualization
summary = [P_ac, Rac, Rdc, ratio];
writematrix(summary, fullfile(pwd, 'problem1_summary.csv'));
disp('Saved problem1_summary.csv for COMSOL post-processing.');

%% Task 4(a): LaTeX formulas printed directly
disp('=== Task 4: LaTeX formulas and numerical solution ===');
disp('delta = sqrt(2/(omega*mu*sigma)) = sqrt(rho/(pi*f*mu))');
disp('J(r) = C*besselj(0,k*r),  k = (1-j)/delta');
disp('C = I*k/(2*pi*a*besselj(1,k*a))');
disp('P_ac = (L/sigma)*2*pi*integral_0^a |J(r)|^2*r dr');
disp('R_ac = P_ac/I^2,  R_dc = L/(sigma*pi*a^2)');

%% Task 4(b): numerical solution with a grid-convergence check
grid_counts = [101, 201, 401, 801, 1601, 3201, 6401];
convergence = zeros(numel(grid_counts), 4);
for q = 1:numel(grid_counts)
    rq = linspace(0, a, grid_counts(q)).';
    Jq = C * besselj(0, k * rq);
    Pq = (L / sigma) * 2 * pi * trapz(rq, abs(Jq).^2 .* rq);
    Raq = Pq / I_rms^2;
    convergence(q, :) = [grid_counts(q), Pq, Raq, Raq / Rdc];
end

writematrix(convergence, fullfile(pwd, 'problem1_grid_convergence.csv'));
disp('Grid convergence columns: Nr, P_ac_W_per_m, R_ac_ohm_per_m, R_ac_over_R_dc');
disp('=== Final numerical solution ===');

% Optional LiveLink sketch (requires COMSOL and LiveLink for MATLAB):
% model = mphopen('problem1_solid_wire.mph');
% model.param.set('f', [num2str(f), '[Hz]']);
% model.param.set('D', [num2str(D), '[m]']);
% model.param.set('I_rms', [num2str(I_rms), '[A]']);
% model.study('std1').run;
% mphsave(model, 'problem1_solid_wire_solved.mph');
