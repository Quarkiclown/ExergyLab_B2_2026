%% Energy and Exergy Analysis: Professional Dashboard (Hot Flow Var)
% Project: Shell-and-Tube Heat Exchanger Performance Analysis

% Author: Al Hadheri, Nistar, Olive
% Date:   April 2026
clear; clc; close all;

%% --- 1. Physical & Geometric Constants ---
% Reference environment and fluid properties
T0  = 283.16;       % Reference temperature for Exergy [K]
k_w = 0.6;          % Thermal conductivity of water [W/mK]
mu_const = 0.0009;  % Dynamic viscosity [Pa.s]
Pr  = 6.966;        % Prandtl Number

% Exchanger Geometry (Shell-and-Tube)
k_s = 15;           % Thermal conductivity of steel [W/mK]
nt  = 6;            % Number of tubes
di  = 0.016;        % Inner tube diameter [m]
do  = 0.018;        % Outer tube diameter [m]
L   = 0.689;        % Heat exchanger length [m]
Ds  = 0.0731;       % Shell diameter [m]

%% --- 2. Input Experimental Data ---
% Counter Flow Experimental Matrix
expDataP = [
% Cols:   Flow_H, Flow_C, Thi,   Tho,   Tci,  Tco
          498.15, 809.82, 28.28, 25.57, 6.72, 8.42;
          792.30, 793.26, 28.48, 26.53, 6.98, 8.84;
          1110.9, 794.34, 28.69, 27.18, 7.52, 9.54;
          1465.5, 800.64, 29.19, 27.96, 8.03, 10.20
        ];

% Counter Flow Experimental Matrix
expDataC = [
% Cols: Flow_H, Flow_C, Thi,   Tho,   Tci,  Tco
        510.75, 808.56, 29.75, 26.82, 7.90, 9.70;
        814.65, 804.24, 29.68, 27.35, 7.76, 10.01;
        1110.9, 807.12, 28.08, 26.41, 7.86, 10.15;
        1449.6, 809.28, 27.93, 26.57, 7.99, 10.40
];

% Extract values into structures for readability
P.flow_h = expDataP(:,1)'; P.flow_c = expDataP(:,2)';
P.Thi = expDataP(:,3)';    P.Tho = expDataP(:,4)';
P.Tci = expDataP(:,5)';    P.Tco = expDataP(:,6)';

C.flow_h = expDataC(:,1)'; C.flow_c = expDataC(:,2)';
C.Thi = expDataC(:,3)';    C.Tho = expDataC(:,4)';
C.Tci = expDataC(:,5)';    C.Tco = expDataC(:,6)';

%% --- 3. Calculations & Modeling ---
[P.I, P.psi, P.epsEx, P.Q] = calcExergy(P.flow_h, P.flow_c, P.Thi, P.Tho, P.Tci, P.Tco, T0);
[C.I, C.psi, C.epsEx, C.Q] = calcExergy(C.flow_h, C.flow_c, C.Thi, C.Tho, C.Tci, C.Tco, T0);

% Simulation Setup: Predictive model across flow range
fh_pred = linspace(min([P.flow_h, C.flow_h])*0.1, 2000, 100); 

% Calculate Mass Flow Ratios (Predicted vs Experimental)
rho_h_mod = getRho(mean([P.Thi, C.Thi])); 
rho_c_mod = getRho(mean([P.Tci, C.Tci]));
mc_const_mod = (mean([P.flow_c, C.flow_c])/3600) * (rho_c_mod/1000);
mh_pred_vec = (fh_pred/3600) * (rho_h_mod/1000);
ratio_pred = mh_pred_vec ./ mc_const_mod;

P_ratio_exp = (P.flow_h .* getRho(P.Thi)) ./ (P.flow_c .* getRho(P.Tci));
C_ratio_exp = (C.flow_h .* getRho(C.Thi)) ./ (C.flow_c .* getRho(C.Tci));

% Execute Heat Exchanger Model (NTU-Effectiveness Method)
[P_mod_I, P_mod_psi, P_mod_epsEx] = runValidatedModel_HotVar(fh_pred, mean(P.flow_c), mean(P.Thi), mean(P.Tci), 'parallel', T0, nt, di, do, L, Ds, k_w, k_s, mu_const, Pr);
[C_mod_I, C_mod_psi, C_mod_epsEx] = runValidatedModel_HotVar(fh_pred, mean(C.flow_c), mean(C.Thi), mean(C.Tci), 'counter', T0, nt, di, do, L, Ds, k_w, k_s, mu_const, Pr);

%% --- 4. Plotting Configuration ---
colorP = [0.12, 0.27, 0.53]; colorC = [0.89, 0.35, 0.13]; opacity = 0.65; 
plotExp = @(ax, x, y, c, m) plot(ax, x, y, '-','Color', c, 'LineWidth', 1.5, 'Marker', m, ...
    'MarkerSize', 8, 'MarkerEdgeColor', c, 'MarkerFaceColor', 'w');
stylePlot = @(ax, yticks, ylims) set(ax, 'LineWidth', 1.1, 'GridAlpha', 0.2, 'TickLabelInterpreter', 'latex', ...
    'FontSize', 11, 'Box', 'on', 'YTick', yticks, 'YLim', ylims, 'XGrid', 'on', 'YGrid', 'on');

%% --- FIGURE 1: Performance vs. Volumetric Hot Flow ---
figure('Name', 'Exergy Analysis: Hot Flow Rate', 'Position', [50, 50, 800, 950], 'Color', 'w');
t1 = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
title(t1, '\textbf{Exergy Performance vs. Hot Flow Rate (Hot Flow Var)}', 'Interpreter', 'latex', 'FontSize', 16);

ax1 = nexttile; hold on;
p_mod1 = plot(fh_pred, P_mod_I, '--', 'Color', [colorP opacity], 'LineWidth', 2);
c_mod1 = plot(fh_pred, C_mod_I, '--', 'Color', [colorC opacity], 'LineWidth', 2);
p_exp1 = plotExp(ax1, P.flow_h, P.I, colorP, 'o'); 
c_exp1 = plotExp(ax1, C.flow_h, C.I, colorC, 's');
ylabel('$\dot{I}$ [W]', 'Interpreter', 'latex'); stylePlot(ax1, 0:100:500, [0 500]);
legend([p_mod1, c_mod1, p_exp1, c_exp1], {'Parallel (Mod)', 'Counter (Mod)', 'Parallel (Exp)', 'Counter (Exp)'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex');

ax2 = nexttile; hold on;
plot(fh_pred, P_mod_psi, '--', 'Color', [colorP opacity], 'LineWidth', 2);
plot(fh_pred, C_mod_psi, '--', 'Color', [colorC opacity], 'LineWidth', 2);
plotExp(ax2, P.flow_h, P.psi, colorP, 'o'); plotExp(ax2, C.flow_h, C.psi, colorC, 's');
ylabel('$\psi$ [\%]', 'Interpreter', 'latex'); stylePlot(ax2, 0:10:50, [0 50]);

ax3 = nexttile; hold on;
plot(fh_pred, P_mod_epsEx, '--', 'Color', [colorP opacity], 'LineWidth', 2);
plot(fh_pred, C_mod_epsEx, '--', 'Color', [colorC opacity], 'LineWidth', 2);
plotExp(ax3, P.flow_h, P.epsEx, colorP, 'o'); plotExp(ax3, C.flow_h, C.epsEx, colorC, 's');
ylabel('$\epsilon_{ex}$ [\%]', 'Interpreter', 'latex'); xlabel('Hot Fluid Flow Rate [l/h]', 'Interpreter', 'latex');
stylePlot(ax3, 0:10:50, [0 50]);
linkaxes([ax1, ax2, ax3], 'x'); xlim(ax1, [0 2000]);

%% --- FIGURE 2: Performance vs. Mass Flow Ratio ---
figure('Name', 'Exergy Analysis: Mass Ratio (Hot Var)', 'Position', [860, 50, 800, 950], 'Color', 'w');
t2 = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
title(t2, '\textbf{Exergy Performance vs. Mass Flow Ratio (Hot Flow Var)}', 'Interpreter', 'latex', 'FontSize', 16);

ax4 = nexttile; hold on;
p_mod2 = plot(ratio_pred, P_mod_I, '--', 'Color', [colorP opacity], 'LineWidth', 2);
c_mod2 = plot(ratio_pred, C_mod_I, '--', 'Color', [colorC opacity], 'LineWidth', 2);
p_exp2 = plotExp(ax4, P_ratio_exp, P.I, colorP, 'o'); 
c_exp2 = plotExp(ax4, C_ratio_exp, C.I, colorC, 's');
ylabel('$\dot{I}$ [W]', 'Interpreter', 'latex'); stylePlot(ax4, 0:100:500, [0 500]);
legend([p_mod2, c_mod2, p_exp2, c_exp2], {'Parallel (Mod)', 'Counter (Mod)', 'Parallel (Exp)', 'Counter (Exp)'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex');

ax5 = nexttile; hold on;
plot(ratio_pred, P_mod_psi, '--', 'Color', [colorP opacity], 'LineWidth', 2);
plot(ratio_pred, C_mod_psi, '--', 'Color', [colorC opacity], 'LineWidth', 2);
plotExp(ax5, P_ratio_exp, P.psi, colorP, 'o'); plotExp(ax5, C_ratio_exp, C.psi, colorC, 's');
ylabel('$\psi$ [\%]', 'Interpreter', 'latex'); stylePlot(ax5, 0:10:50, [0 50]);

ax6 = nexttile; hold on;
plot(ratio_pred, P_mod_epsEx, '--', 'Color', [colorP opacity], 'LineWidth', 2);
plot(ratio_pred, C_mod_epsEx, '--', 'Color', [colorC opacity], 'LineWidth', 2);
plotExp(ax6, P_ratio_exp, P.epsEx, colorP, 'o'); plotExp(ax6, C_ratio_exp, C.epsEx, colorC, 's');
ylabel('$\epsilon_{ex}$ [\%]', 'Interpreter', 'latex'); xlabel('Mass Flow Ratio $\dot{m}_h / \dot{m}_c$', 'Interpreter', 'latex');
stylePlot(ax6, 0:10:50, [0 50]);
linkaxes([ax4, ax5, ax6], 'x'); xlim(ax4, [0 2.5]);

%% --- HELPER FUNCTIONS ---
function rho = getRho(T)
    % Density correlation for water [kg/m^3]
    rho = -0.0000001*T.^4 + 0.00004*T.^3 - 0.0074*T.^2 + 0.0477*T + 999.91;
end
function cp = getCp(T)
    % Specific heat correlation for water [J/kgK]
    cp = (0.000000003*T.^4 - 0.0000007*T.^3 + 0.00007*T.^2 - 0.0028*T + 4.2173)*1000;
end

function [I, psi, epsilon_ex, Q] = calcExergy(fh, fc, Thi, Tho, Tci, Tco, T0)
    % Thermodynamic Exergy Balance
    rho_h = getRho((Thi+Tho)/2); cp_h = getCp((Thi+Tho)/2);
    rho_c = getRho((Tci+Tco)/2); cp_c = getCp((Tci+Tco)/2);
    mh = (fh/3600).*(rho_h/1000); mc = (fc/3600).*(rho_c/1000);
    ThiK = Thi+273.15; ThoK = Tho+273.15; TciK = Tci+273.15; TcoK = Tco+273.15;
    
    Q = mh.*cp_h.*(ThiK - ThoK); % Heat transfer [W]
    dExH = mh.*cp_h.*((ThiK - ThoK) - T0.*log(ThiK./ThoK)); % Exergy supplied
    dExC = mc.*cp_c.*((TcoK - TciK) - T0.*log(TcoK./TciK)); % Exergy recovered
    
    I = dExH - dExC;                % Exergy Destruction [W]
    psi = abs((dExC ./ dExH)) * 100; % Rational Exergetic Efficiency [%]
    
    Ch = mh.*cp_h; Cc = mc.*cp_c; Cmin = min(Ch, Cc);
    epsilon_ex = (Ch.*(ThiK - ThoK)) ./ (Cmin.*(ThiK - TciK)) * 100; % Thermal Effectiveness [%]
    
    % Data Cleaning
    I(isnan(I)) = 0; psi(isnan(psi)) = 0; epsilon_ex(isnan(epsilon_ex)) = 0;
end

function [I, psi, epsilon_ex] = runValidatedModel_HotVar(fh_list, fc, Thi, Tci, mode, T0, nt, di, do, L, Ds, k_w, k_s, mu, Pr)
    % Predictive Model using Convection Correlations and NTU method
    I = zeros(size(fh_list)); psi = zeros(size(fh_list)); epsilon_ex = zeros(size(fh_list));
    Area_shell = (pi/4)*(Ds^2 - nt*do^2); Perimeter_shell = pi*Ds + nt*pi*do; Dh = 4*Area_shell/Perimeter_shell;
    Ai = nt*pi*di*L; Ao = nt*pi*do*L; R_wall = log(do/di)/(2*pi*k_s*L*nt);
    
    for i = 1:length(fh_list)
        rho_h = getRho(Thi); cp_h = getCp(Thi); rho_c = getRho(Tci); cp_c = getCp(Tci);
        mh = (fh_list(i)/3600)*(rho_h/1000); mc = (fc/3600)*(rho_c/1000);
        
        % Nusselt & Reynolds Correlations (Physics Engine)
        hi = (0.023*((rho_c*(mc/(rho_c*nt*pi/4*di^2))*di)/mu)^0.8 * Pr^0.4 * (1+6*di/L) * k_w)/di;
        ho = (0.36*((rho_h*(mh/(rho_h*Area_shell))*Dh)/mu)^0.55 * Pr^(1/3) * k_w)/Dh;
        UA = 1/(1/(hi*Ai) + R_wall + 1/(ho*Ao));
        
        % NTU-Effectiveness Branching Logic
        Cc = mc*cp_c; Ch = mh*cp_h; Cmin = min(Cc,Ch); Cr = Cmin/max(Cc,Ch); NTU = UA/Cmin;
        if strcmp(mode,'parallel')
            eps = (1-exp(-NTU*(1+Cr)))/(1+Cr);
        else
            if Cr==1, eps = NTU/(1+NTU); 
            else, eps = (1-exp(-NTU*(1-Cr)))/(1-Cr*exp(-NTU*(1-Cr))); end
        end
        
        % Predict Outlets and calculate 2nd Law Metrics
        Q = eps*Cmin*((Thi+273.15)-(Tci+273.15));
        ThoK = (Thi+273.15)-Q/Ch; TcoK = (Tci+273.15)+Q/Cc;
        [I(i), psi(i), epsilon_ex(i), ~] = calcExergy(fh_list(i), fc, Thi, ThoK-273.15, Tci, TcoK-273.15, T0);
    end
end
