%% Energy and Exergy Analysis: Professional Dashboard (Variable Properties)
% Scenario: Constant Cold Flow, Varying Hot Flow (0 - 2000 l/h)
clear; clc; close all;

%% --- 1. Geometric & Physical Constants ---
T0  = 283.16; k_w = 0.6; 
mu_const = 0.001; % Dynamic viscosity used for Reynolds
Pr  = 6.966; k_s = 16; nt  = 6; di  = 0.016; 
do  = 0.018; L   = 0.689; Ds  = 0.0731;

%% --- 2. Input Experimental Data (Varying Hot Flow) ---
P.flow_h = [498.15, 792.30, 1110.90, 1465.50];
P.flow_c = [809.82, 793.26, 794.34, 800.64];
P.Thi = [28.28, 28.48, 28.69, 29.19];    P.Tho = [25.57, 26.53, 27.18, 27.96];
P.Tci = [6.72, 6.98, 7.52, 8.03];        P.Tco = [8.42, 8.84, 9.54, 10.20];

C.flow_h = [510.75, 814.65, 1110.90, 1449.60];
C.flow_c = [808.56, 804.24, 807.12, 809.28];
C.Thi = [29.75, 29.68, 28.08, 27.93];    C.Tho = [26.82, 27.35, 26.41, 26.57];
C.Tci = [7.90, 7.76, 7.86, 7.99];        C.Tco = [9.70, 10.01, 10.15, 10.40];

%% --- 3. Calculations ---
% Added P.Q and C.Q outputs to avoid table errors later
[P.I, P.psi, P.epsEx, P.Q] = calcExergy(P.flow_h, P.flow_c, P.Thi, P.Tho, P.Tci, P.Tco, T0);
[C.I, C.psi, C.epsEx, C.Q] = calcExergy(C.flow_h, C.flow_c, C.Thi, C.Tho, C.Tci, C.Tco, T0);

fh_pred = linspace(10, 2000, 100); 
[P_mod_I, P_mod_psi, P_mod_epsEx] = runValidatedModel(fh_pred, mean(P.flow_c), ...
    mean(P.Thi), mean(P.Tci), 'parallel', T0, nt, di, do, L, Ds, k_w, k_s, mu_const, Pr);
[C_mod_I, C_mod_psi, C_mod_epsEx] = runValidatedModel(fh_pred, mean(C.flow_c), ...
    mean(C.Thi), mean(C.Tci), 'counter', T0, nt, di, do, L, Ds, k_w, k_s, mu_const, Pr);

% Ratio: Hot / Cold (Since Hot is the varying parameter)
P_ratio_exp = P.flow_h ./ P.flow_c;
C_ratio_exp = C.flow_h ./ C.flow_c;
ratio_pred  = fh_pred ./ mean([P.flow_c, C.flow_c]);

%% --- 4. Plotting & Visuals ---
colorP = [0.12, 0.27, 0.53]; % Midnight Blue
colorC = [0.89, 0.35, 0.13]; % Sunset Orange
bgGray = [0.99, 0.99, 0.99]; 
stylePlot = @(ax) set(ax, 'LineWidth', 1.1, 'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 11, 'Box', 'on', 'Color', bgGray);
expPlot = @(x, y, c, m) plot(x, y, 'Color', c, 'LineWidth', 1.5, ...
    'Marker', m, 'MarkerSize', 8, 'MarkerEdgeColor', c, 'MarkerFaceColor', 'w');

%% --- FIGURE 1: Exergy vs Hot Flow Rate ---
figure('Name', 'Exergy Analysis: Hot Flow Rate Variation', 'Position', [50, 50, 800, 950], 'Color', 'w');
t1 = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

ax1_1 = nexttile; hold on; grid on;
p1 = plot(fh_pred, P_mod_I, '--', 'Color', [colorP, 0.65], 'LineWidth', 1.2); 
p2 = plot(fh_pred, C_mod_I, '--', 'Color', [colorC, 0.65], 'LineWidth', 1.2);
e1 = expPlot(P.flow_h, P.I, colorP, 'o');
e2 = expPlot(C.flow_h, C.I, colorC, 's');
ylabel('$\mathbf{\dot{I}}$ \textbf{[W]}', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Exergy Performance vs. Hot Flow Rate}', 'Interpreter', 'latex', 'FontSize', 15);
stylePlot(ax1_1); xticklabels(ax1_1, {}); 
lgd1 = legend([p1, p2, e1, e2], {'Parallel (Model)', 'Counter (Model)', 'Parallel (Exp)', 'Counter (Exp)'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 10);
lgd1.Layout.Tile = 'north'; 

ax1_2 = nexttile; hold on; grid on;
plot(fh_pred, P_mod_psi, '--', 'Color', [colorP, 0.65], 'LineWidth', 1.2);
plot(fh_pred, C_mod_psi, '--', 'Color', [colorC, 0.65], 'LineWidth', 1.2);
expPlot(P.flow_h, P.psi, colorP, 'o');
expPlot(C.flow_h, C.psi, colorC, 's');
ylabel('$\mathbf{\psi}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
stylePlot(ax1_2); xticklabels(ax1_2, {}); 

ax1_3 = nexttile; hold on; grid on;
plot(fh_pred, P_mod_epsEx, '--', 'Color', [colorP, 0.65], 'LineWidth', 1.2);
plot(fh_pred, C_mod_epsEx, '--', 'Color', [colorC, 0.65], 'LineWidth', 1.2);
expPlot(P.flow_h, P.epsEx, colorP, 'o');
expPlot(C.flow_h, C.epsEx, colorC, 's');
ylabel('$\mathbf{\epsilon_{ex}}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('\textbf{Hot Fluid Flow Rate, } $\mathbf{\dot{V}_h}$ \textbf{[l/h]}', 'Interpreter', 'latex', 'FontSize', 12);
stylePlot(ax1_3);
linkaxes([ax1_1, ax1_2, ax1_3], 'x'); xlim(ax1_1, [0 2000]);

%% --- FIGURE 2: Exergy vs Flow Ratio (mh/mc) ---
figure('Name', 'Exergy Analysis: Qhot Var (Flow Ratio)', 'Position', [860, 50, 800, 950], 'Color', 'w');
t2 = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

ax2_1 = nexttile; hold on; grid on;
p3 = plot(ratio_pred, P_mod_I, '--', 'Color', [colorP, 0.65], 'LineWidth', 1.2);
p4 = plot(ratio_pred, C_mod_I, '--', 'Color', [colorC, 0.65], 'LineWidth', 1.2);
[~, sP] = sort(P_ratio_exp); [~, sC] = sort(C_ratio_exp);
e3 = expPlot(P_ratio_exp(sP), P.I(sP), colorP, 'o');
e4 = expPlot(C_ratio_exp(sC), C.I(sC), colorC, 's');
ylabel('$\mathbf{\dot{I}}$ \textbf{[W]}', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Exergy Performance vs. Mass Flow Ratio ($\dot{m}_h/\dot{m}_c$)}', 'Interpreter', 'latex', 'FontSize', 15);
stylePlot(ax2_1); xticklabels(ax2_1, {}); 
lgd2 = legend([p3, p4, e3, e4], {'Parallel (Model)', 'Counter (Model)', 'Parallel (Exp)', 'Counter (Exp)'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 10);
lgd2.Layout.Tile = 'north'; 

ax2_2 = nexttile; hold on; grid on;
plot(ratio_pred, P_mod_psi, '--', 'Color', [colorP, 0.65], 'LineWidth', 1.2);
plot(ratio_pred, C_mod_psi, '--', 'Color', [colorC, 0.65], 'LineWidth', 1.2);
expPlot(P_ratio_exp(sP), P.psi(sP), colorP, 'o');
expPlot(C_ratio_exp(sC), C.psi(sC), colorC, 's');
ylabel('$\mathbf{\psi}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
stylePlot(ax2_2); xticklabels(ax2_2, {}); 

ax2_3 = nexttile; hold on; grid on;
plot(ratio_pred, P_mod_epsEx, '--', 'Color', [colorP, 0.65], 'LineWidth', 1.2);
plot(ratio_pred, C_mod_epsEx, '--', 'Color', [colorC, 0.65], 'LineWidth', 1.2);
expPlot(P_ratio_exp(sP), P.epsEx(sP), colorP, 'o');
expPlot(C_ratio_exp(sC), C.epsEx(sC), colorC, 's');
ylabel('$\mathbf{\epsilon_{ex}}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('\textbf{Mass Flow Rate Ratio, } $\mathbf{\dot{m}_h / \dot{m}_c}$', 'Interpreter', 'latex', 'FontSize', 12);
stylePlot(ax2_3);
linkaxes([ax2_1, ax2_2, ax2_3], 'x'); xlim(ax2_1, [0 2.5]);

%% --- Display Combined Results Table ---
fprintf('\n==========================================================\n');
fprintf('       EXPERIMENTAL RESULTS: ENERGY & EXERGY SUMMARY\n');
fprintf('==========================================================\n');
VarNames = {'HotFlow_lh', 'ColdFlow_lh', 'Heat_Q_W', 'ExDest_I_W', 'Efficiency_pct'};
T_P = table(P.flow_h', P.flow_c', P.Q', P.I', P.psi', 'VariableNames', VarNames);
fprintf('\n[Parallel Flow Configuration]:\n');
disp(T_P);
T_C = table(C.flow_h', C.flow_c', C.Q', C.I', C.psi', 'VariableNames', VarNames);
fprintf('\n[Counter Flow Configuration]:\n');
disp(T_C);

%% --- HELPER FUNCTIONS ---
% Density and Cp Calculations at each point
function rho = getRho(T_avg)
    rho = -0.0000001*(T_avg).^4 + 0.00004*(T_avg).^3 - 0.0074*(T_avg).^2 + 0.0477*(T_avg) + 999.91;
end

function cp = getCp(T_avg)
    cp = (0.000000003*(T_avg).^4 - 0.0000007*(T_avg).^3 + 0.00007*(T_avg).^2 - 0.0028*(T_avg) + 4.2173)*1000;
end

% Numerical simulation of NTU method mathmatical model 
function [I, psi, epsilon_ex, Q] = calcExergy(fh, fc, Thi, Tho, Tci, Tco, T0)
    Tavg_h = (Thi + Tho) / 2;
    Tavg_c = (Tci + Tco) / 2;
    rho_h = getRho(Tavg_h); cp_h = getCp(Tavg_h);
    rho_c = getRho(Tavg_c); cp_c = getCp(Tavg_c);
    
    mh = (fh / 3600) .* (rho_h / 1000); mc = (fc / 3600) .* (rho_c / 1000);
    ThiK = Thi+273.15; ThoK = Tho+273.15; TciK = Tci+273.15; TcoK = Tco+273.15;
    
    % Heat transfer calculation
    Q = mh .* cp_h .* (ThiK - ThoK); 
    
    dExH = mh .* cp_h .* ((ThiK - ThoK) - T0 .* log(ThiK./ThoK)); 
    dExC = mc .* cp_c .* ((TcoK - TciK) - T0 .* log(TcoK./TciK)); 
    dEx_max = min(mh, mc) .* cp_h .* ((ThiK - TciK) - T0 .* log(ThiK./TciK));
    
    I = dExH - dExC; 
    psi = abs((dExC ./ dExH)) * 100;
    epsilon_ex = (dExH ./ dEx_max) * 100;
    psi(isnan(psi)) = 0; epsilon_ex(isnan(epsilon_ex)) = 0;
end

function [I, psi, epsilon_ex] = runValidatedModel(fh_list, fc, Thi, Tci, mode, T0, nt, di, do, L, Ds, k_w, k_s, mu, Pr)
    I = zeros(size(fh_list)); psi = zeros(size(fh_list)); epsilon_ex = zeros(size(fh_list));
    ThiK = Thi + 273.15; TciK = Tci + 273.15;
    
    for i = 1:length(fh_list)
        rho_h = getRho(Thi); cp_h = getCp(Thi);
        rho_c = getRho(Tci); cp_c = getCp(Tci);
        
        mh = (fh_list(i)/3600)*(rho_h/1000);
        mc = (fc/3600)*(rho_c/1000);
        
        Dh = (Ds^2 - nt*do^2) / (Ds + nt*do);
        u_s = mh / (rho_h * (pi/4 * (Ds^2 - nt*do^2))); Re_s = (u_s * Dh * rho_h) / mu;
        Nu_s = 0.36 * Re_s^0.55 * Pr^(1/3); ho = (Nu_s * k_w) / Dh;
        
        u_t = mc / (rho_c * nt * (pi/4 * di^2)); Re_t = (u_t * di * rho_c) / mu; 
        Nu_t = 0.023 * Re_t^0.8 * Pr^0.4 * (1 + 6*di/L); hi = (Nu_t * k_w) / di;
        
        R_wall = log(do/di) / (2 * pi * k_s * L * nt);
        US = 1 / (1/(ho * nt*pi*do*L) + R_wall + 1/(hi * nt*pi*di*L));
        
        Cc = mc*cp_c; Ch = mh*cp_h; Cmin = min(Cc,Ch); Cr = Cmin/max(Cc,Ch);
        NTU = US / Cmin;
        
        if strcmp(mode, 'parallel')
            eps_hx = (1 - exp(-NTU*(1+Cr))) / (1+Cr);
        else
            eps_hx = (1 - exp(-NTU*(1-Cr))) / (1 - Cr*exp(-NTU*(1-Cr)));
        end
        
        Q = eps_hx * Cmin * (ThiK - TciK);
        ThoK = ThiK - Q/Ch; TcoK = TciK + Q/Cc;
        
        [I(i), psi(i), epsilon_ex(i)] = calcExergy(fh_list(i), fc, Thi, ThoK-273.15, Tci, TcoK-273.15, T0);
    end
end