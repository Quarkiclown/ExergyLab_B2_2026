%% Compressor_Exergy_Lab_AlHadheri_Nistar_Olive.m
% This script performs analysis of compressor experimental data and
% improved theoretical modeling for energy and exergy evaluation.
%
% It defines thermodynamic constants, loads two experimental cases, computes
% compressor performance metrics, fits isentropic efficiency models, and
% evaluates energy/exergy performance metrics.
%
% Author: Al Hadheri, Nistar, Olive (Year 4, Semester 8 Exergy Analysis Lab)
% Date: April 3, 2026

clear; clc; close all;

%% 1. Thermodynamic constants
R_kJ   = 0.287;
cp_kJ  = 1.004;
gamma  = 1.4;
Cd     = 0.61;
S1     = 0.000581;
S2     = 0.000123;
T0     = 298.15;

R  = R_kJ * 1000;
cp = cp_kJ * 1000;
beta = sqrt(S2 / S1);
D_pipe = sqrt(4 * S1 / pi);

%% 2. Experimental data

% CASE 1
data1.delta_P1 = [34, 32, 47, 47, 47, 46];
data1.P2_bar   = [1, 2, 3, 4, 5, 6];
data1.P4_mbar  = [992, 992, 993, 993, 993, 993];
data1.T1_C     = [27.2, 27.1, 27.4, 27.5, 27.6, 28.0];
data1.T2_C     = [51.8, 53.0, 54.0, 54.7, 56.3, 58.1];
data1.T3_C     = [22.3, 22.3, 22.4, 22.4, 22.5, 22.7];
data1.W_act    = [105, 146, 172, 196, 219, 245];
data1.RPM      = 750 * ones(size(data1.P2_bar));

% CASE 2
data2.delta_P1 = [44, 43, 44, 28];
data2.P2_bar   = [3, 3, 3, 3];
data2.P4_mbar  = [993, 993, 994, 994];
data2.T1_C     = [28.3, 28.9, 28.9, 28.8];
data2.T2_C     = [60.0, 61.8, 66.2, 75.4];
data2.T3_C     = [22.9, 22.9, 23.0, 23.4];
data2.W_act    = [118, 154, 196, 242];
data2.RPM      = [500, 700, 900, 1100];

%% 3. Process experimental data
r1 = analyseExperimentalCase(data1, R, cp, gamma, Cd, S1, S2, beta, D_pipe, T0);
r2 = analyseExperimentalCase(data2, R, cp, gamma, Cd, S1, S2, beta, D_pipe, T0);

%% 4. FIT η_is (energy efficiency) & apply to theory

% ---- CASE 1 ----
eta_is_exp1 = r1.eta_ise / 100;   % fraction
PR1 = r1.PR;

p1 = polyfit(PR1, eta_is_exp1, 2);
eta_is_fit1 = polyval(p1, PR1);

% ---- CASE 2 ----
eta_is_exp2 = r2.eta_ise / 100;
RPM2 = data2.RPM;

p2 = polyfit(RPM2, eta_is_exp2, 2);
eta_is_fit2 = polyval(p2, RPM2);

%% 5. Compute theoretical efficiencies for BOTH energy and exergy
r1_th = theoreticalFromFittedEfficiencyExergy(r1, eta_is_fit1);
r2_th = theoreticalFromFittedEfficiencyExergy(r2, eta_is_fit2);

%% 6. Plot experimental vs theoretical

% ---- Case 1 ----
figure('Color','w');
plot(data1.P2_bar, r1.eta_energy, 'o-', 'LineWidth',2); hold on;
plot(data1.P2_bar, r1.eta_exergy, 's-', 'LineWidth',2);
plot(data1.P2_bar, r1_th.eta_energy, '--', 'LineWidth',2);
plot(data1.P2_bar, r1_th.eta_exergy, '--', 'LineWidth',2);
grid on;
xlabel('P2 [bar]');
ylabel('Efficiency [%]');
title('Case 1: Variable η_{is} (energy + exergy)');
legend('η_{energy} exp','η_{exergy} exp','η_{energy} theory','η_{exergy} theory');

% ---- Case 2 ----
figure('Color','w');
plot(data2.RPM, r2.eta_energy, 'o-', 'LineWidth',2); hold on;
plot(data2.RPM, r2.eta_exergy, 's-', 'LineWidth',2);
plot(data2.RPM, r2_th.eta_energy, '--', 'LineWidth',2);
plot(data2.RPM, r2_th.eta_exergy, '--', 'LineWidth',2);
grid on;
xlabel('RPM');
ylabel('Efficiency [%]');
title('Case 2: Variable η_{is} (energy + exergy)');
legend('η_{energy} exp','η_{exergy} exp','η_{energy} theory','η_{exergy} theory');

%% =============================================================================
function r = analyseExperimentalCase(d, R, cp, gamma, Cd, S1, S2, beta, D_pipe, T0)
n = numel(d.delta_P1);

r.mdot = zeros(1,n);
r.Q_vol = zeros(1,n);
r.PR = zeros(1,n);
r.W_ise = zeros(1,n);
r.W_iso = zeros(1,n);
r.eta_ise = zeros(1,n);
r.eta_iso = zeros(1,n);
r.eta_energy = zeros(1,n);
r.eta_exergy = zeros(1,n);

for k = 1:n
    T1 = d.T1_C(k)+273.15;
    T2 = d.T2_C(k)+273.15;
    T3 = d.T3_C(k)+273.15;

    P_in = d.P4_mbar(k)*100;
    P_out = d.P2_bar(k)*1e5 + P_in;

    rho = P_in/(R*T3);

    Qv = Cd*S2*sqrt(2*d.delta_P1(k)/rho)/sqrt(1-beta^4);
    mdot = rho*Qv;

    PR = P_out/P_in;

    W_ise = mdot*cp*T1*(PR^((gamma-1)/gamma)-1);
    W_iso = mdot*R*T1*log(PR);

    eta_ise = (W_ise/d.W_act(k))*100;
    eta_iso = (W_iso/d.W_act(k))*100;

    r.Q_vol(k) = Qv;
    r.mdot(k) = mdot;
    r.PR(k) = PR;
    r.W_ise(k) = W_ise;
    r.W_iso(k) = W_iso;
    r.eta_ise(k) = eta_ise;
    r.eta_iso(k) = eta_iso;
    r.eta_energy(k) = eta_ise;
    r.eta_exergy(k) = eta_iso;
end
end

%% =============================================================================
function rth = theoreticalFromFittedEfficiencyExergy(r_exp, eta_is_fit)
% Computes both energy and exergy efficiencies for theory
rth.W_act_th   = r_exp.W_ise ./ eta_is_fit;
rth.eta_energy = eta_is_fit * 100;                       % energy
rth.eta_exergy = (r_exp.W_iso ./ rth.W_act_th) * 100;    % exergy
end
