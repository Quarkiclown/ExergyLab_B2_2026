clear; clc; close all;
% ---------------------------------------------------------
% Compressor Energy and Exergy Analysis - Dual Case Study
% ---------------------------------------------------------

%% 1. Thermodynamic Constants (Shared)
R = 0.287;          % Specific gas constant for air [kJ/kg-K]
cp = 1.004;         % Specific heat at constant pressure [kJ/kg-K]
gamma = 1.4;        % Ratio of specific heats
Cd = 0.61;          % Discharge coefficient
S1 = 0.000581;      % Pipe cross section [m^2]
S2 = 0.000123;      % Orifice cross section [m^2]
T0 = 298.15;        % Dead state temperature (25 C) [K]

%% 2. Experimental Data Input
% --- CASE 1: Fixed Motor RPM (750), Varying Pressure ---
data1.delta_P1 = [34, 32, 47, 47, 47, 46];        
data1.P2_bar   = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];    
data1.P4_mbar  = [992, 992, 993, 993, 993, 993];   
data1.T1_C     = [27.2, 27.1, 27.4, 27.5, 27.6, 28.0];
data1.T2_C     = [51.8, 53.0, 54.0, 54.7, 56.3, 58.1];
data1.T3_C     = [22.3, 22.3, 22.4, 22.4, 22.5, 22.7];
data1.W_act    = [105, 146, 172, 196, 219, 245];     

% --- CASE 2: Fixed Pressure (3.0 bar), Varying Frequency ---
data2.delta_P1 = [44, 43, 44, 28];        
data2.P2_bar   = [3.0, 3.0, 3.0, 3.0]; % Fixed pressure    
data2.P4_mbar  = [993, 993, 994, 994];   
data2.T1_C     = [28.3, 28.9, 28.9, 28.8];
data2.T2_C     = [60.0, 61.8, 66.2, 75.4];
data2.T3_C     = [22.9, 22.9, 23.0, 23.4];
data2.W_act    = [118, 154, 196, 242];   
% IMPORTANT: Add your actual recorded RPM values here!
data2.RPM      = [500, 700, 900, 1100]; 

%% 3. Mathematical Processing
all_data = {data1, data2};
case_names = {'Case 1: Fixed RPM', 'Case 2: Fixed Pressure'};

% Create empty cell arrays to store results for plotting later
m_dot_all = cell(1, 2);
eta_en_all = cell(1, 2);
eta_ex_all = cell(1, 2);
X_dest_exp_all = cell(1, 2);
X_dest_theo_all = cell(1, 2); % New array for theoretical exergy destruction

%% 3. Mathematical Processing
all_data = {data1, data2};
case_names = {'Case 1: Fixed RPM', 'Case 2: Fixed Pressure'};

m_dot_all = cell(1, 2);
eta_en_all = cell(1, 2);
eta_ex_all = cell(1, 2);
X_dest_exp_all = cell(1, 2);
X_dest_theo_all = cell(1, 2);

for i = 1:2
    d = all_data{i};
    
    % Unit Conversions
    T1 = d.T1_C + 273.15; 
    T2 = d.T2_C + 273.15;
    T3 = d.T3_C + 273.15;
    P_atm = d.P4_mbar * 100;
    P_in  = P_atm;
    P_out = (d.P2_bar * 10^5) + P_atm;
    P0_orifice = P_atm + d.delta_P1; 
    
    % Flow Rate Calculations
    rho0 = P0_orifice ./ ((R * 1000) .* T3); 
    S_ratio_term = (S1^2 / S2^2) - 1;
    Q = Cd * S1 * sqrt((2 * d.delta_P1) ./ (rho0 * S_ratio_term)); 
    m_dot_all{i} = rho0 .* Q; 
    
    % Energy Analysis
    W_ise = m_dot_all{i} * cp * 1000 .* T1 .* ((P_out./P_in).^((gamma-1)/gamma) - 1);
    eta_en_all{i} = (W_ise ./ d.W_act) * 100;
    
    % Exergy Analysis
    term1 = cp * (T2 - T1);
    term2 = T0 * (cp * log(T2./T1) - R * log(P_out./P_in));
    delta_psi = m_dot_all{i} .* (term1 - term2) * 1000; 
    
    % 1. Experimental Exergy Destruction
    X_dest_exp_all{i} = d.W_act - delta_psi;
    eta_ex_all{i} = (delta_psi ./ d.W_act) * 100;
    
    % 2. Theoretical Exergy Destruction (Gouy-Stodola) — FIXED
    % Isentropic outlet temperature (ideal, reversible reference)
    T2s = T1 .* (P_out ./ P_in) .^ ((gamma - 1) / gamma);
    
    % Entropy generation = excess heating above isentropic path (always >= 0)
    S_gen = cp * log(T2 ./ T2s);  % kJ/kg-K
    
    % Gouy-Stodola: X_dest = T0 * m_dot * S_gen
    X_dest_theo_all{i} = T0 * m_dot_all{i} .* S_gen * 1000;  % Watts
end

%% 4. FIGURE 1: Efficiency Analysis (2x2 Grid)
figure(1);
set(gcf, 'Name', 'Efficiency Analysis', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
tlo1 = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'loose');
title(tlo1, 'Thermodynamic Efficiency Analysis', 'FontSize', 16, 'FontWeight', 'bold');

% Tile 1: Case 1 - Efficiency vs Pressure
nexttile;
plot(data1.P2_bar, eta_en_all{1}, '-o', 'LineWidth', 2); hold on;
plot(data1.P2_bar, eta_ex_all{1}, '-s', 'LineWidth', 2);
title('Case 1: Efficiency vs Pressure Ratio');
xlabel('$P_2$ (bar)', 'Interpreter', 'latex'); ylabel('Efficiency (\%)', 'Interpreter', 'latex');
grid on; legend('$\eta_{en}$', '$\eta_{ex}$', 'Interpreter', 'latex', 'Location', 'best');

% Tile 2: Case 1 - Efficiency vs Mass Flow Rate
nexttile;
plot(m_dot_all{1}*1e4, eta_en_all{1}, '-o', 'LineWidth', 2); hold on;
plot(m_dot_all{1}*1e4, eta_ex_all{1}, '-s', 'LineWidth', 2);
title('Case 1: Efficiency vs Mass Flow Rate');
xlabel('$\dot{m} \times 10^{-4}$ (kg/s)', 'Interpreter', 'latex'); ylabel('Efficiency (\%)', 'Interpreter', 'latex');
grid on; legend('$\eta_{en}$', '$\eta_{ex}$', 'Interpreter', 'latex', 'Location', 'best');

% Tile 3: Case 2 - Efficiency vs RPM
nexttile;
plot(data2.RPM, eta_en_all{2}, '-^', 'Color', [0.46 0.67 0.18], 'LineWidth', 2); hold on;
plot(data2.RPM, eta_ex_all{2}, '-d', 'Color', [0.49 0.18 0.55], 'LineWidth', 2);
title('Case 2: Efficiency vs Motor RPM');
xlabel('Motor Speed (RPM)'); ylabel('Efficiency (\%)', 'Interpreter', 'latex');
grid on; legend('$\eta_{en}$', '$\eta_{ex}$', 'Interpreter', 'latex', 'Location', 'best');

% Tile 4: Case 2 - Efficiency vs Mass Flow Rate
nexttile;
plot(m_dot_all{2}*1e4, eta_en_all{2}, '-^', 'Color', [0.46 0.67 0.18], 'LineWidth', 2); hold on;
plot(m_dot_all{2}*1e4, eta_ex_all{2}, '-d', 'Color', [0.49 0.18 0.55], 'LineWidth', 2);
title('Case 2: Efficiency vs Mass Flow Rate');
xlabel('$\dot{m} \times 10^{-4}$ (kg/s)', 'Interpreter', 'latex'); ylabel('Efficiency (\%)', 'Interpreter', 'latex');
grid on; legend('$\eta_{en}$', '$\eta_{ex}$', 'Interpreter', 'latex', 'Location', 'best');

%% 5. FIGURE 2: Exergy Destruction Analysis (2x2 Grid) SUPERIMPOSED
figure(2);
set(gcf, 'Name', 'Exergy Destruction', 'Units', 'normalized', 'Position', [0.15 0.15 0.8 0.8]);
tlo2 = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'loose');
title(tlo2, 'Experimental vs Theoretical Exergy Destruction', 'FontSize', 16, 'FontWeight', 'bold');

% Tile 1: Case 1 - X_dest vs Pressure
nexttile;
plot(data1.P2_bar, X_dest_exp_all{1}, '-o', 'Color', '#D95319', 'LineWidth', 2, 'MarkerFaceColor', 'w'); hold on;
plot(data1.P2_bar, X_dest_theo_all{1}, '--k', 'LineWidth', 1.5, 'Marker', '*');
title('Case 1: \chi_{dest} vs Pressure Ratio');
xlabel('$P_2$ (bar)', 'Interpreter', 'latex'); ylabel('$\dot{X}_{dest}$ (W)', 'Interpreter', 'latex');
grid on; legend('Experimental (Total)', 'Theoretical (Internal)', 'Location', 'best');

% Tile 2: Case 1 - X_dest vs Mass Flow Rate
nexttile;
plot(m_dot_all{1}*1e4, X_dest_exp_all{1}, '-o', 'Color', '#D95319', 'LineWidth', 2, 'MarkerFaceColor', 'w'); hold on;
plot(m_dot_all{1}*1e4, X_dest_theo_all{1}, '--k', 'LineWidth', 1.5, 'Marker', '*');
title('Case 1: \chi_{dest} vs Mass Flow Rate');
xlabel('$\dot{m} \times 10^{-4}$ (kg/s)', 'Interpreter', 'latex'); ylabel('$\dot{X}_{dest}$ (W)', 'Interpreter', 'latex');
grid on; legend('Experimental (Total)', 'Theoretical (Internal)', 'Location', 'best');

% Tile 3: Case 2 - X_dest vs RPM
nexttile;
plot(data2.RPM, X_dest_exp_all{2}, '-^', 'Color', '#0072BD', 'LineWidth', 2, 'MarkerFaceColor', 'w'); hold on;
plot(data2.RPM, X_dest_theo_all{2}, '--k', 'LineWidth', 1.5, 'Marker', '*');
title('Case 2: \chi_{dest} vs Motor RPM');
xlabel('Motor Speed (RPM)'); ylabel('$\dot{X}_{dest}$ (W)', 'Interpreter', 'latex');
grid on; legend('Experimental (Total)', 'Theoretical (Internal)', 'Location', 'best');

% Tile 4: Case 2 - X_dest vs Mass Flow Rate
nexttile;
plot(m_dot_all{2}*1e4, X_dest_exp_all{2}, '-^', 'Color', '#0072BD', 'LineWidth', 2, 'MarkerFaceColor', 'w'); hold on;
plot(m_dot_all{2}*1e4, X_dest_theo_all{2}, '--k', 'LineWidth', 1.5, 'Marker', '*');
title('Case 2: \chi_{dest} vs Mass Flow Rate');
xlabel('$\dot{m} \times 10^{-4}$ (kg/s)', 'Interpreter', 'latex'); ylabel('$\dot{X}_{dest}$ (W)', 'Interpreter', 'latex');
grid on; legend('Experimental (Total)', 'Theoretical (Internal)', 'Location', 'best');