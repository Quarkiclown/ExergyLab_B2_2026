%% Energy and Exergy Analysis: Professional Dashboard (Variable Properties)
% Scenario: Constant Hot Flow, Varying Cold Flow (0 - 2000 l/h)
clear; clc; close all;

%% --- 1. Physical & Geometric Constants ---
% Static cp and rho removed to use temperature-dependent equations
T0  = 283.16; k_w = 0.6; 
mu_const = 0.001; % Dynamic viscosity for Reynolds calc
Pr  = 6.966; k_s = 16; nt  = 6; di  = 0.016; 
do  = 0.018; L   = 0.689; Ds  = 0.0731;

%% --- 2. Input Experimental Data ---
P.flow_h = [823.05, 838.2, 838.95, 831]; 
P.flow_c = [447.83, 757.97, 1099.25, 1386];
P.Thi = [28.96, 29.10, 29.42, 29.18];   P.Tho = [27.47, 27.23, 26.98, 26.88];
P.Tci = [6.42, 6.76, 7, 7.29];          P.Tco = [9.08, 8.73, 8.73, 8.7];

C.flow_h = [803.25, 800.25, 811.65, 811.95]; 
C.flow_c = [480.96, 777.78, 1099.44, 1398.42];
C.Thi = [28.51, 29.17, 29.14, 29.5];    C.Tho = [26.94, 26.75, 26.41, 26.39];
C.Tci = [7.06, 7.06, 7.18, 7.76];       C.Tco = [9.74, 9.5, 9.21, 9.53];

%% --- 3. Calculations ---
[P.I, P.psi, P.epsEx, P.Q] = calcExergy(P.flow_h, P.flow_c, P.Thi, P.Tho, P.Tci, P.Tco, T0);
[C.I, C.psi, C.epsEx, C.Q] = calcExergy(C.flow_h, C.flow_c, C.Thi, C.Tho, C.Tci, C.Tco, T0);

fc_pred = linspace(10, 2000, 100); 
[P_mod_I, P_mod_psi, P_mod_epsEx, P_mod_Q] = runValidatedModel(fc_pred, mean(P.flow_h), ...
    mean(P.Thi), mean(P.Tci), 'parallel', T0, nt, di, do, L, Ds, k_w, k_s, mu_const, Pr);
[C_mod_I, C_mod_psi, C_mod_epsEx, C_mod_Q] = runValidatedModel(fc_pred, mean(C.flow_h), ...
    mean(C.Thi), mean(C.Tci), 'counter', T0, nt, di, do, L, Ds, k_w, k_s, mu_const, Pr);

P_ratio_exp = P.flow_h ./ P.flow_c;
C_ratio_exp = C.flow_h ./ C.flow_c;
ratio_pred  = mean([P.flow_h, C.flow_h]) ./ fc_pred;

%% --- 4. Plotting & Visuals ---
colorP = [0.12, 0.27, 0.53]; % Midnight Blue
colorC = [0.89, 0.35, 0.13]; % Sunset Orange
bgGray = [0.99, 0.99, 0.99]; 
stylePlot = @(ax) set(ax, 'LineWidth', 1.1, 'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 11, 'Box', 'on', 'Color', bgGray);
mScatter = @(x, y, c, m) scatter(x, y, 70, m, 'LineWidth', 1.8, ...
    'MarkerEdgeColor', c, 'MarkerFaceColor', 'w');

%% --- FIGURE 1: Original Dashboard (vs Flow Rate) ---
figure('Name', 'Exergy Analysis: Flow Rate', 'Position', [50, 50, 800, 950], 'Color', 'w');
t1 = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
ax1_1 = nexttile; hold on; grid on;
p1 = plot(fc_pred, P_mod_I, '--', 'Color', [colorP, 0.65], 'LineWidth', 2); 
p2 = plot(fc_pred, C_mod_I, '--', 'Color', [colorC, 0.65], 'LineWidth', 2); 
plot(P.flow_c, P.I, '-', 'Color', colorP, 'LineWidth', 1.5); 
plot(C.flow_c, C.I, '-', 'Color', colorC, 'LineWidth', 1.5);
e1 = mScatter(P.flow_c, P.I, colorP, 'o');
e2 = mScatter(C.flow_c, C.I, colorC, 's');
ylabel('$\mathbf{\dot{I}}$ \textbf{[W]}', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Exergy Performance vs. Cold Flow Rate}', 'Interpreter', 'latex', 'FontSize', 15);
stylePlot(ax1_1); xticklabels(ax1_1, {}); 
lgd1 = legend([p1, p2, e1, e2], {'Parallel (Model)', 'Counter (Model)', 'Parallel (Exp)', 'Counter (Exp)'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 10);
lgd1.Layout.Tile = 'north'; 

ax1_2 = nexttile; hold on; grid on;
plot(fc_pred, P_mod_psi, '--', 'Color', [colorP, 0.65], 'LineWidth', 2); 
plot(fc_pred, C_mod_psi, '--', 'Color', [colorC, 0.65], 'LineWidth', 2);
plot(P.flow_c, P.psi, '-', 'Color', colorP, 'LineWidth', 1.5);
plot(C.flow_c, C.psi, '-', 'Color', colorC, 'LineWidth', 1.5);
mScatter(P.flow_c, P.psi, colorP, 'o');
mScatter(C.flow_c, C.psi, colorC, 's');
ylabel('$\mathbf{\psi}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
stylePlot(ax1_2); xticklabels(ax1_2, {}); 

ax1_3 = nexttile; hold on; grid on;
plot(fc_pred, P_mod_epsEx, '--', 'Color', [colorP, 0.65], 'LineWidth', 2);
plot(fc_pred, C_mod_epsEx, '--', 'Color', [colorC, 0.65], 'LineWidth', 2);
plot(P.flow_c, P.epsEx, '-', 'Color', colorP, 'LineWidth', 1.5);
plot(C.flow_c, C.epsEx, '-', 'Color', colorC, 'LineWidth', 1.5);
mScatter(P.flow_c, P.epsEx, colorP, 'o');
mScatter(C.flow_c, C.epsEx, colorC, 's');
ylabel('$\mathbf{\epsilon_{ex}}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('\textbf{Cold Fluid Flow Rate, } $\mathbf{\dot{V}_c}$ \textbf{[l/h]}', 'Interpreter', 'latex', 'FontSize', 12);
stylePlot(ax1_3);
linkaxes([ax1_1, ax1_2, ax1_3], 'x'); xlim(ax1_1, [0 2000]);

%% --- FIGURE 2: New Dashboard (vs Flow Ratio mh/mc) ---
figure('Name', 'Exergy Analysis: Qcold Var (Flow Ratio)', 'Position', [860, 50, 800, 950], 'Color', 'w');
t2 = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
ax2_1 = nexttile; hold on; grid on;
p3 = plot(ratio_pred, P_mod_I, '--', 'Color', [colorP, 0.65], 'LineWidth', 2); 
p4 = plot(ratio_pred, C_mod_I, '--', 'Color', [colorC, 0.65], 'LineWidth', 2); 
plot(P_ratio_exp, P.I, '-', 'Color', colorP, 'LineWidth', 1.5); 
plot(C_ratio_exp, C.I, '-', 'Color', colorC, 'LineWidth', 1.5);
e3 = mScatter(P_ratio_exp, P.I, colorP, 'o');
e4 = mScatter(C_ratio_exp, C.I, colorC, 's');
ylabel('$\mathbf{\dot{I}}$ \textbf{[W]}', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Exergy Performance vs. Mass Flow Ratio}', 'Interpreter', 'latex', 'FontSize', 15);
stylePlot(ax2_1); xticklabels(ax2_1, {}); 
lgd2 = legend([p3, p4, e3, e4], {'Parallel (Model)', 'Counter (Model)', 'Parallel (Exp)', 'Counter (Exp)'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 10);
lgd2.Layout.Tile = 'north'; 

ax2_2 = nexttile; hold on; grid on;
plot(ratio_pred, P_mod_psi, '--', 'Color', [colorP, 0.65], 'LineWidth', 2);
plot(ratio_pred, C_mod_psi, '--', 'Color', [colorC, 0.65], 'LineWidth', 2);
plot(P_ratio_exp, P.psi, '-', 'Color', colorP, 'LineWidth', 1.5);
plot(C_ratio_exp, C.psi, '-', 'Color', colorC, 'LineWidth', 1.5);
mScatter(P_ratio_exp, P.psi, colorP, 'o');
mScatter(C_ratio_exp, C.psi, colorC, 's');
ylabel('$\mathbf{\psi}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
stylePlot(ax2_2); xticklabels(ax2_2, {}); 

ax2_3 = nexttile; hold on; grid on;
plot(ratio_pred, P_mod_epsEx, '--', 'Color', [colorP, 0.65], 'LineWidth', 2);
plot(ratio_pred, C_mod_epsEx, '--', 'Color', [colorC, 0.65], 'LineWidth', 2);
plot(P_ratio_exp, P.epsEx, '-', 'Color', colorP, 'LineWidth', 1.5);
plot(C_ratio_exp, C.epsEx, '-', 'Color', colorC, 'LineWidth', 1.5);
mScatter(P_ratio_exp, P.epsEx, colorP, 'o');
mScatter(C_ratio_exp, C.epsEx, colorC, 's');
ylabel('$\mathbf{\epsilon_{ex}}$ \textbf{[\%]}', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('\textbf{Flow Rate Ratio, } $\mathbf{\dot{m}_h / \dot{m}_c}$', 'Interpreter', 'latex', 'FontSize', 12);
stylePlot(ax2_3);
linkaxes([ax2_1, ax2_2, ax2_3], 'x'); xlim(ax2_1, [0 4.5]); 

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

% Experemental Data Calculations
function [I, psi, epsilon_ex, Q] = calcExergy(fh, fc, Thi, Tho, Tci, Tco, T0)
    Tavg_h = (Thi + Tho) / 2;
    Tavg_c = (Tci + Tco) / 2;
    rho_h = getRho(Tavg_h); cp_h = getCp(Tavg_h);
    rho_c = getRho(Tavg_c); cp_c = getCp(Tavg_c);
    mh = (fh / 3600) .* (rho_h / 1000); 
    mc = (fc / 3600) .* (rho_c / 1000);
    ThiK = Thi+273.15; ThoK = Tho+273.15; TciK = Tci+273.15; TcoK = Tco+273.15;
    Q = mh .* cp_h .* (ThiK - ThoK); 
    dExH = mh .* cp_h .* ((ThiK - ThoK) - T0 .* log(ThiK./ThoK)); 
    dExC = mc .* cp_c .* ((TcoK - TciK) - T0 .* log(TcoK./TciK)); 
    dEx_max = min(mh, mc) .* cp_h .* ((ThiK - TciK) - T0 .* log(ThiK./TciK));
    I = dExH - dExC; 
    psi = abs((dExC ./ dExH)) * 100;
    epsilon_ex = (dExH ./ dEx_max) * 100;
    psi(isnan(psi)) = 0; epsilon_ex(isnan(epsilon_ex)) = 0;
end

% Numerical simulation of NTU method mathmatical model 
function [I, psi, epsilon_ex, Q_list] = runValidatedModel(fc_list, fh, Thi, Tci, mode, T0, nt, di, do, L, Ds, k_w, k_s, mu, Pr)
    I = zeros(size(fc_list)); psi = zeros(size(fc_list)); 
    epsilon_ex = zeros(size(fc_list)); Q_list = zeros(size(fc_list));
    ThiK = Thi + 273.15; TciK = Tci + 273.15;
    for i = 1:length(fc_list)
        rho_h = getRho(Thi); cp_h = getCp(Thi);
        rho_c = getRho(Tci); cp_c = getCp(Tci);
        mh = (fh/3600)*(rho_h/1000);
        mc = (fc_list(i)/3600)*(rho_c/1000);
        u_t = mc / (rho_c * nt * (pi/4 * di^2)); Re_t = (u_t * di * rho_c) / mu; 
        Nu_t = 0.023 * Re_t^0.8 * Pr^0.4 * (1 + 6*di/L); hi = (Nu_t * k_w) / di;
        Dh = (Ds^2 - nt*do^2) / (Ds + nt*do);
        u_s = mh / (rho_h * (pi/4 * (Ds^2 - nt*do^2))); Re_s = (u_s * Dh * rho_h) / mu;
        Nu_s = 0.36 * Re_s^0.55 * Pr^(1/3); ho = (Nu_s * k_w) / Dh;
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
        Q_list(i) = Q;
        ThoK = ThiK - Q/Ch; TcoK = TciK + Q/Cc;
        [I(i), psi(i), epsilon_ex(i), ~] = calcExergy(fh, fc_list(i), Thi, ThoK-273.15, Tci, TcoK-273.15, T0);
    end
end