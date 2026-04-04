%% =========================================================================
%  EXPERIMENTAL vs THEORETICAL COMPARISON
%  Shell-and-tube Water/Water Heat Exchanger
%  Metrics: Qh, İ (exergy destruction), Sgen, ε (effectiveness), ψ (exergy efficiency)
%  4 datasets: counter/parallel × cold-varying/hot-varying
%
%  Theory: continuous curves swept from 0 to 2000 l/h
%  Experiment: discrete points overlaid on the theoretical curves
%
%  ψ = İ / ΔẊh × 100 [%] — fraction of hot-side exergy destroyed.
%  Valid for all conditions including sub-ambient cold inlet (Tc,in < T0).
% =========================================================================
clear; clc; close all;

%% =========================================================================
%  SECTION 1 — FLUID PROPERTY FUNCTIONS  (polynomial fits, T in °C)
% =========================================================================
rho_w = @(T) -1e-7*T.^4 + 4e-5*T.^3 - 7.4e-3*T.^2 + 4.77e-2*T + 999.91;
cp_w  = @(T) (3e-9*T.^4 - 7e-7*T.^3 + 7e-5*T.^2 - 2.8e-3*T + 4.2173)*1000;
mu_w  = @(T) (3e-5*T.^4 - 8.5e-3*T.^3 + 9.355e-1*T.^2 - 53.679*T + 1765.4)*1e-6;
lam_w = @(T) (6e-7*T.^4 - 1e-4*T.^3 - 1.1e-3*T.^2 + 1.9509*T + 560.85)/1000;
Pr_w  = @(T) mu_w(T).*cp_w(T)./lam_w(T);

%% =========================================================================
%  SECTION 2 — GEOMETRY & CONSTANTS
% =========================================================================
di  = 0.016;   de  = 0.018;   nt  = 6;
L   = 0.689;   Ds  = 0.0731;  k_w = 15;
T0  = 283.16;  % dead-state temperature [K]

Sm   = pi*di^2/4;
Se   = pi*Ds^2/4;
A_sh = Se - nt*Sm;
Dh   = 4*A_sh / (pi*Ds + nt*pi*de);
Si   = nt*pi*di*L;
Se_a = nt*pi*de*L;

%% =========================================================================
%  SECTION 3 — EXPERIMENTAL DATASETS
%  Columns: [Qvh(l/h), Qvc(l/h), Thin(°C), Thout(°C), Tcin(°C), Tcout(°C)]
% =========================================================================

% Counter-flow, cold flow varying (~800 l/h hot fixed)
data_CC = [
    803.25,  480.96,  28.51,  26.94,  7.06,   9.74;
    800.25,  777.78,  29.17,  26.75,  7.06,   9.50;
    811.65, 1099.44,  29.14,  26.41,  7.18,   9.21;
    811.95, 1398.42,  29.50,  26.39,  7.76,   9.53;
];

% Parallel-flow, cold flow varying
data_PC = [
    823.05,  447.84,  28.96,  27.47,  6.42,   9.08;
    838.20,  757.98,  29.10,  27.23,  6.76,   8.73;
    838.95, 1099.26,  29.42,  26.98,  7.00,   8.73;
    831.00, 1386.00,  29.18,  26.88,  7.29,   8.70;
];

% Counter-flow, hot flow varying (~808 l/h cold fixed)
data_CH = [
     510.75,  808.56,  29.75,  26.82,  7.90,   9.70;
     814.65,  804.24,  29.68,  27.35,  7.76,  10.01;
    1110.90,  807.12,  28.08,  26.41,  7.86,  10.15;
    1449.60,  809.28,  27.93,  26.57,  7.99,  10.40;
];

% Parallel-flow, hot flow varying
data_PH = [
     498.15,  809.82,  28.28,  25.57,  6.72,   8.42;
     792.30,  793.26,  28.48,  26.53,  6.98,   8.84;
    1110.90,  794.34,  28.69,  27.18,  7.52,   9.54;
    1465.50,  800.64,  29.19,  27.96,  8.03,  10.20;
];

%% =========================================================================
%  SECTION 4 — NOMINAL INLET CONDITIONS FOR THEORETICAL SWEEPS
%  Use the mean of the measured inlet temperatures across each dataset
% =========================================================================
Thin_cold_C  = mean(data_CC(:,3));   % mean hot inlet, cold-varying datasets
Tcin_cold_C  = mean(data_CC(:,5));   % mean cold inlet, cold-varying datasets
Thin_cold_P  = mean(data_PC(:,3));
Tcin_cold_P  = mean(data_PC(:,5));

Thin_hot_C   = mean(data_CH(:,3));   % mean hot inlet, hot-varying datasets
Tcin_hot_C   = mean(data_CH(:,5));
Thin_hot_P   = mean(data_PH(:,3));
Tcin_hot_P   = mean(data_PH(:,5));

Qvh_fixed_C  = mean(data_CC(:,1));   % fixed hot flow for cold-sweep (counter)
Qvh_fixed_P  = mean(data_PC(:,1));   % fixed hot flow for cold-sweep (parallel)
Qvc_fixed_C  = mean(data_CH(:,2));   % fixed cold flow for hot-sweep (counter)
Qvc_fixed_P  = mean(data_PH(:,2));   % fixed cold flow for hot-sweep (parallel)

%% =========================================================================
%  SECTION 5 — COMPUTE EXPERIMENTAL VALUES FROM MEASURED TEMPERATURES
% =========================================================================

function res = compute_exp(data, T0, rho_w, cp_w)
    n = size(data, 1);
    Qh    = zeros(n,1);
    Idest = zeros(n,1);
    Sgen  = zeros(n,1);
    eps   = zeros(n,1);
    psi   = zeros(n,1);

    for i = 1:n
        Qvh_lh=data(i,1); Qvc_lh=data(i,2);
        Thin=data(i,3);   Thout=data(i,4);
        Tcin=data(i,5);   Tcout=data(i,6);

        Qvh=Qvh_lh/3.6e6;  Qvc=Qvc_lh/3.6e6;
        Tfh=(Thin+Thout)/2; Tfc=(Tcin+Tcout)/2;

        mh=Qvh*rho_w(Tfh); mc=Qvc*rho_w(Tfc);
        cph=cp_w(Tfh);      cpc=cp_w(Tfc);
        Ch=mh*cph;          Cc=mc*cpc;
        Cmin=min(Ch,Cc);

        ThinK=Thin+273.15; ThoutK=Thout+273.15;
        TcinK=Tcin+273.15; TcoutK=Tcout+273.15;

        Qh(i)    = mh*cph*(ThinK-ThoutK);
        dXh      = mh*(cph*(ThinK-ThoutK) - T0*cph*log(ThinK/ThoutK));
        dXc      = mc*(cpc*(TcoutK-TcinK) - T0*cpc*log(TcoutK/TcinK));
        Idest(i) = dXh - dXc;
        Sgen(i)  = mh*cph*log(ThinK/ThoutK) + mc*cpc*log(TcoutK/TcinK);
        eps(i)   = Ch*(ThinK-ThoutK)/(Cmin*(ThinK-TcinK))*100;
        psi(i)   = Idest(i)/dXh*100;
    end
    res.Qh    = Qh;
    res.Idest = Idest;
    res.Sgen  = Sgen*1000;   % → mW/K
    res.eps   = eps;
    res.psi   = psi;
end

%% =========================================================================
%  SECTION 6 — NTU-EFFECTIVENESS THEORETICAL MODEL (single operating point)
% =========================================================================

function res = compute_theory_point(Qvh_lh, Qvc_lh, Thin_C, Tcin_C, config, T0, ...
    rho_w, cp_w, mu_w, lam_w, Pr_w, di, de, nt, L, k_w, Sm, A_sh, Dh, Si, Se_a)

    Qvh=Qvh_lh/3.6e6;  Qvc=Qvc_lh/3.6e6;
    ThinK=Thin_C+273.15; TcinK=Tcin_C+273.15;
    ThoutK=ThinK-2;      TcoutK=TcinK+2;

    for iter = 1:80
        Tfh=(ThinK+ThoutK)/2-273.15; Tfc=(TcinK+TcoutK)/2-273.15;
        mh=Qvh*rho_w(Tfh);  mc=Qvc*rho_w(Tfc);
        cph=cp_w(Tfh);  cpc=cp_w(Tfc);
        muh=mu_w(Tfh);  muc=mu_w(Tfc);
        lh=lam_w(Tfh);  lc=lam_w(Tfc);
        Ch=mh*cph;  Cc=mc*cpc;
        Cmin=min(Ch,Cc);  Cr=Cmin/max(Ch,Cc);

        % Shell side (hot): Nu = 0.36 Re^0.55 Pr^(1/3)
        V_sh=Qvh/A_sh;
        Re_h=V_sh*Dh/(muh/rho_w(Tfh));
        Nu_h=0.36*Re_h^0.55*Pr_w(Tfh)^(1/3);
        he=Nu_h*lh/Dh;

        % Tube side (cold): Nu = 0.023 Re^0.8 Pr^0.4 (1+6D/L)
        V_tb=Qvc/(nt*Sm);
        Re_c=V_tb*di/(muc/rho_w(Tfc));
        Nu_c=0.023*Re_c^0.8*Pr_w(Tfc)^0.4*(1+6*di/L);
        hi=Nu_c*lc/di;

        Req=1/(he*Se_a)+log(de/di)/(2*pi*k_w*L)+1/(hi*Si);
        UA=1/Req;  NTU=UA/Cmin;

        if strcmp(config,'counter')
            if abs(Cr-1)<1e-6
                eff=NTU/(1+NTU);
            else
                eff=(1-exp(-NTU*(1-Cr)))/(1-Cr*exp(-NTU*(1-Cr)));
            end
        else
            eff=(1-exp(-NTU*(1+Cr)))/(1+Cr);
        end

        Q=eff*Cmin*(ThinK-TcinK);
        Th2=ThinK-Q/Ch;  Tc2=TcinK+Q/Cc;
        if abs(Th2-ThoutK)<1e-5 && abs(Tc2-TcoutK)<1e-5
            break;
        end
        ThoutK=Th2;  TcoutK=Tc2;
    end

    ThC=ThoutK-273.15;  TcC=TcoutK-273.15;
    mh=Qvh*rho_w((Thin_C+ThC)/2);   mc=Qvc*rho_w((Tcin_C+TcC)/2);
    cph=cp_w((Thin_C+ThC)/2);        cpc=cp_w((Tcin_C+TcC)/2);
    Ch=mh*cph;  Cc=mc*cpc;  Cmin=min(Ch,Cc);

    Qh    = mh*cph*(ThinK-ThoutK);
    dXh   = mh*(cph*(ThinK-ThoutK)-T0*cph*log(ThinK/ThoutK));
    dXc   = mc*(cpc*(TcoutK-TcinK)-T0*cpc*log(TcoutK/TcinK));
    Idest = dXh-dXc;
    Sgen  = mh*cph*log(ThinK/ThoutK)+mc*cpc*log(TcoutK/TcinK);
    eps_v = Ch*(ThinK-ThoutK)/(Cmin*(ThinK-TcinK))*100;
    psi_v = Idest/dXh*100;

    res.Qh    = Qh;
    res.Idest = Idest;
    res.Sgen  = Sgen*1000;
    res.eps   = eps_v;
    res.psi   = psi_v;
end

%% =========================================================================
%  SECTION 7 — SWEEP THEORETICAL MODEL OVER 0–2000 l/h
% =========================================================================

Q_sweep = linspace(10, 2000, 300);   % avoid exactly 0 (div-by-zero in Re)
n_sw    = length(Q_sweep);

% Pre-allocate: fields for counter and parallel, cold-sweep and hot-sweep
fields = {'Qh','Idest','Sgen','eps','psi'};
for f = fields
    Th_cc_cold.(f{1}) = zeros(n_sw,1);
    Th_pc_cold.(f{1}) = zeros(n_sw,1);
    Th_cc_hot.(f{1})  = zeros(n_sw,1);
    Th_ph_hot.(f{1})  = zeros(n_sw,1);
end

fprintf('Computing theoretical curves (0–2000 l/h)...\n');
for k = 1:n_sw
    q = Q_sweep(k);

    % Cold-flow sweep: vary Qvc, fix Qvh at mean measured value
    r = compute_theory_point(Qvh_fixed_C, q, Thin_cold_C, Tcin_cold_C, 'counter', T0, ...
        rho_w,cp_w,mu_w,lam_w,Pr_w, di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    for f = fields, Th_cc_cold.(f{1})(k) = r.(f{1}); end

    r = compute_theory_point(Qvh_fixed_P, q, Thin_cold_P, Tcin_cold_P, 'parallel', T0, ...
        rho_w,cp_w,mu_w,lam_w,Pr_w, di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    for f = fields, Th_pc_cold.(f{1})(k) = r.(f{1}); end

    % Hot-flow sweep: vary Qvh, fix Qvc at mean measured value
    r = compute_theory_point(q, Qvc_fixed_C, Thin_hot_C, Tcin_hot_C, 'counter', T0, ...
        rho_w,cp_w,mu_w,lam_w,Pr_w, di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    for f = fields, Th_cc_hot.(f{1})(k) = r.(f{1}); end

    r = compute_theory_point(q, Qvc_fixed_P, Thin_hot_P, Tcin_hot_P, 'parallel', T0, ...
        rho_w,cp_w,mu_w,lam_w,Pr_w, di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    for f = fields, Th_ph_hot.(f{1})(k) = r.(f{1}); end
end
fprintf('Done.\n\n');

%% =========================================================================
%  SECTION 8 — COMPUTE EXPERIMENTAL VALUES
% =========================================================================

E_CC = compute_exp(data_CC, T0, rho_w, cp_w);
E_PC = compute_exp(data_PC, T0, rho_w, cp_w);
E_CH = compute_exp(data_CH, T0, rho_w, cp_w);
E_PH = compute_exp(data_PH, T0, rho_w, cp_w);

% X-axis positions of the experimental points
x_cold_C = data_CC(:,2);
x_cold_P = data_PC(:,2);
x_hot_C  = data_CH(:,1);
x_hot_P  = data_PH(:,1);

%% =========================================================================
%  SECTION 9 — PLOT HELPERS
% =========================================================================
blue  = [0.216 0.541 0.867];
coral = [0.847 0.353 0.188];
lw_th = 2.0;    % theory line width
lw_ex = 1.5;    % experimental line (connecting markers)
ms    = 9;      % marker size

% Theory = solid continuous line
% Experiment = markers only (no connecting line between discrete points)
function plot_theory_exp(Q_sw, y_th_C, y_th_P, x_e_C, y_e_C, x_e_P, y_e_P, ...
                          blue, coral, lw_th, ms, ylbl)
    % Theory curves
    plot(Q_sw, y_th_C, '-',  'Color',blue,  'LineWidth',lw_th, ...
         'DisplayName','Counter (theory)');
    plot(Q_sw, y_th_P, '--', 'Color',coral, 'LineWidth',lw_th, ...
         'DisplayName','Parallel (theory)');
    % Experimental scatter
    scatter(x_e_C, y_e_C, 70, blue,  'o', 'filled', ...
            'DisplayName','Counter (exp.)', 'LineWidth',1);
    scatter(x_e_P, y_e_P, 70, coral, 's', 'filled', ...
            'DisplayName','Parallel (exp.)', 'LineWidth',1);
    grid on; box on;
    ylabel(ylbl, 'FontSize',11);
    set(gca,'FontSize',10);
    xlim([0 2000]);
end

%% =========================================================================
%  FIGURE 1 — COLD FLOW VARYING (0–2000 l/h)
% =========================================================================
metric_ylabels = {
    'Qh (W)';
    'Exergy destruction \dot{I} (W)';
    '\dot{S}_{gen} (mW/K)';
    'Effectiveness \epsilon (%)';
    '\psi = \dot{I}/\Delta\dot{X}_h (%)';
};

th_cold_C_vals = {Th_cc_cold.Qh, Th_cc_cold.Idest, Th_cc_cold.Sgen, Th_cc_cold.eps, Th_cc_cold.psi};
th_cold_P_vals = {Th_pc_cold.Qh, Th_pc_cold.Idest, Th_pc_cold.Sgen, Th_pc_cold.eps, Th_pc_cold.psi};
ex_cold_C_vals = {E_CC.Qh, E_CC.Idest, E_CC.Sgen, E_CC.eps, E_CC.psi};
ex_cold_P_vals = {E_PC.Qh, E_PC.Idest, E_PC.Sgen, E_PC.eps, E_PC.psi};

figure('Name','Cold flow varying — theory 0–2000 l/h', ...
       'NumberTitle','off','Position',[50 50 1350 920]);
sgtitle({'Cold flow rate varying  |  Hot fixed at mean measured value', ...
         sprintf('Q_{v,hot} = %.0f l/h (counter),  %.0f l/h (parallel)', ...
                 Qvh_fixed_C, Qvh_fixed_P), ...
         sprintf('T_{h,in} = %.1f°C,  T_{c,in} = %.1f°C (counter nominal)', ...
                 Thin_cold_C, Tcin_cold_C)}, ...
        'FontSize',12,'FontWeight','bold');

for k = 1:5
    subplot(2,3,k);
    hold on;
    plot_theory_exp(Q_sweep, th_cold_C_vals{k}, th_cold_P_vals{k}, ...
                    x_cold_C, ex_cold_C_vals{k}, x_cold_P, ex_cold_P_vals{k}, ...
                    blue, coral, lw_th, ms, metric_ylabels{k});
    xlabel('Cold flow rate (l/h)', 'FontSize',11);
    title(metric_ylabels{k}, 'FontSize',11);
    if k == 1
        legend('Location','best','FontSize',9);
    end
end

subplot(2,3,6);
axis off;
annotation_text = {
    '\bf Legend', ...
    '\color[rgb]{0.216 0.541 0.867}— Counter-flow (theory)', ...
    '\color[rgb]{0.847 0.353 0.188}– – Parallel-flow (theory)', ...
    '\color[rgb]{0.216 0.541 0.867}\circ Counter-flow (experimental)', ...
    '\color[rgb]{0.847 0.353 0.188}\square Parallel-flow (experimental)', ...
    '', ...
    'Theory: NTU-\epsilon method, continuous 0–2000 l/h', ...
    'Experiment: 4 discrete operating points per config.'
};
text(0.05, 0.95, annotation_text, 'Units','normalized', ...
     'VerticalAlignment','top','FontSize',10,'Interpreter','tex');

%% =========================================================================
%  FIGURE 2 — HOT FLOW VARYING (0–2000 l/h)
% =========================================================================
th_hot_C_vals = {Th_cc_hot.Qh, Th_cc_hot.Idest, Th_cc_hot.Sgen, Th_cc_hot.eps, Th_cc_hot.psi};
th_hot_P_vals = {Th_ph_hot.Qh, Th_ph_hot.Idest, Th_ph_hot.Sgen, Th_ph_hot.eps, Th_ph_hot.psi};
ex_hot_C_vals = {E_CH.Qh, E_CH.Idest, E_CH.Sgen, E_CH.eps, E_CH.psi};
ex_hot_P_vals = {E_PH.Qh, E_PH.Idest, E_PH.Sgen, E_PH.eps, E_PH.psi};

figure('Name','Hot flow varying — theory 0–2000 l/h', ...
       'NumberTitle','off','Position',[100 80 1350 920]);
sgtitle({'Hot flow rate varying  |  Cold fixed at mean measured value', ...
         sprintf('Q_{v,cold} = %.0f l/h (counter),  %.0f l/h (parallel)', ...
                 Qvc_fixed_C, Qvc_fixed_P), ...
         sprintf('T_{h,in} = %.1f°C,  T_{c,in} = %.1f°C (counter nominal)', ...
                 Thin_hot_C, Tcin_hot_C)}, ...
        'FontSize',12,'FontWeight','bold');

for k = 1:5
    subplot(2,3,k);
    hold on;
    plot_theory_exp(Q_sweep, th_hot_C_vals{k}, th_hot_P_vals{k}, ...
                    x_hot_C, ex_hot_C_vals{k}, x_hot_P, ex_hot_P_vals{k}, ...
                    blue, coral, lw_th, ms, metric_ylabels{k});
    xlabel('Hot flow rate (l/h)', 'FontSize',11);
    title(metric_ylabels{k}, 'FontSize',11);
    if k == 1
        legend('Location','best','FontSize',9);
    end
end

subplot(2,3,6);
axis off;
text(0.05, 0.95, annotation_text, 'Units','normalized', ...
     'VerticalAlignment','top','FontSize',10,'Interpreter','tex');

%% =========================================================================
%  FIGURE 3 — PARITY PLOTS (exp vs theory at matched operating points)
% =========================================================================

% Compute theory at the exact experimental flow rates for parity comparison
n_CC=size(data_CC,1); n_PC=size(data_PC,1);
n_CH=size(data_CH,1); n_PH=size(data_PH,1);

T_CC_pts = struct('Qh',zeros(n_CC,1),'Idest',zeros(n_CC,1));
T_PC_pts = struct('Qh',zeros(n_PC,1),'Idest',zeros(n_PC,1));
T_CH_pts = struct('Qh',zeros(n_CH,1),'Idest',zeros(n_CH,1));
T_PH_pts = struct('Qh',zeros(n_PH,1),'Idest',zeros(n_PH,1));

for i=1:n_CC
    r=compute_theory_point(data_CC(i,1),data_CC(i,2),data_CC(i,3),data_CC(i,5),'counter',T0,...
        rho_w,cp_w,mu_w,lam_w,Pr_w,di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    T_CC_pts.Qh(i)=r.Qh; T_CC_pts.Idest(i)=r.Idest;
end
for i=1:n_PC
    r=compute_theory_point(data_PC(i,1),data_PC(i,2),data_PC(i,3),data_PC(i,5),'parallel',T0,...
        rho_w,cp_w,mu_w,lam_w,Pr_w,di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    T_PC_pts.Qh(i)=r.Qh; T_PC_pts.Idest(i)=r.Idest;
end
for i=1:n_CH
    r=compute_theory_point(data_CH(i,1),data_CH(i,2),data_CH(i,3),data_CH(i,5),'counter',T0,...
        rho_w,cp_w,mu_w,lam_w,Pr_w,di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    T_CH_pts.Qh(i)=r.Qh; T_CH_pts.Idest(i)=r.Idest;
end
for i=1:n_PH
    r=compute_theory_point(data_PH(i,1),data_PH(i,2),data_PH(i,3),data_PH(i,5),'parallel',T0,...
        rho_w,cp_w,mu_w,lam_w,Pr_w,di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
    T_PH_pts.Qh(i)=r.Qh; T_PH_pts.Idest(i)=r.Idest;
end

figure('Name','Parity plots','NumberTitle','off','Position',[150 150 950 420]);
sgtitle('Parity: experimental vs theoretical (at matched operating points)', ...
        'FontSize',13,'FontWeight','bold');

% İ parity
subplot(1,2,1); hold on; grid on; box on;
all_Ie = [E_CC.Idest; E_PC.Idest; E_CH.Idest; E_PH.Idest];
all_It = [T_CC_pts.Idest; T_PC_pts.Idest; T_CH_pts.Idest; T_PH_pts.Idest];
scatter(E_CC.Idest, T_CC_pts.Idest, 70, blue,  'o','filled','DisplayName','Counter cold-var');
scatter(E_PC.Idest, T_PC_pts.Idest, 70, coral, 'o','filled','DisplayName','Parallel cold-var');
scatter(E_CH.Idest, T_CH_pts.Idest, 70, blue,  's','filled','DisplayName','Counter hot-var');
scatter(E_PH.Idest, T_PH_pts.Idest, 70, coral, 's','filled','DisplayName','Parallel hot-var');
lims=[min([all_Ie;all_It])*0.93, max([all_Ie;all_It])*1.05];
plot(lims,lims,'k--','LineWidth',1.5,'DisplayName','Perfect agreement');
xlim(lims); ylim(lims);
xlabel('İ experimental (W)','FontSize',12); ylabel('İ theoretical (W)','FontSize',12);
title('Exergy destruction İ','FontSize',12);
legend('Location','southeast','FontSize',9); axis square;

% Qh parity
subplot(1,2,2); hold on; grid on; box on;
all_Qe = [E_CC.Qh; E_PC.Qh; E_CH.Qh; E_PH.Qh];
all_Qt = [T_CC_pts.Qh; T_PC_pts.Qh; T_CH_pts.Qh; T_PH_pts.Qh];
scatter(E_CC.Qh, T_CC_pts.Qh, 70, blue,  'o','filled','DisplayName','Counter cold-var');
scatter(E_PC.Qh, T_PC_pts.Qh, 70, coral, 'o','filled','DisplayName','Parallel cold-var');
scatter(E_CH.Qh, T_CH_pts.Qh, 70, blue,  's','filled','DisplayName','Counter hot-var');
scatter(E_PH.Qh, T_PH_pts.Qh, 70, coral, 's','filled','DisplayName','Parallel hot-var');
lims=[min([all_Qe;all_Qt])*0.90, max([all_Qe;all_Qt])*1.07];
plot(lims,lims,'k--','LineWidth',1.5,'DisplayName','Perfect agreement');
xlim(lims); ylim(lims);
xlabel('Qh experimental (W)','FontSize',12); ylabel('Qh theoretical (W)','FontSize',12);
title('Heat transfer Qh','FontSize',12);
legend('Location','southeast','FontSize',9); axis square;

%% =========================================================================
%  FIGURE 4 — EXERGY EFFICIENCY ψ ONLY (0–2000 l/h, large readable plot)
% =========================================================================
figure('Name','Exergy efficiency psi','NumberTitle','off','Position',[200 200 1050 450]);
sgtitle({'\psi = \dot{I}/\Delta\dot{X}_h \times 100  [% of hot-side exergy destroyed]', ...
         'Lower \psi = less irreversibility = better thermodynamic performance'}, ...
        'FontSize',12,'FontWeight','bold');

subplot(1,2,1); hold on;
plot(Q_sweep, Th_cc_cold.psi, '-',  'Color',blue,  'LineWidth',lw_th, ...
     'DisplayName','Counter (theory)');
plot(Q_sweep, Th_pc_cold.psi, '--', 'Color',coral, 'LineWidth',lw_th, ...
     'DisplayName','Parallel (theory)');
scatter(x_cold_C, E_CC.psi, 80, blue,  'o','filled','DisplayName','Counter (exp.)');
scatter(x_cold_P, E_PC.psi, 80, coral, 's','filled','DisplayName','Parallel (exp.)');
xlabel('Cold flow rate (l/h)','FontSize',12);
ylabel('\psi (%)','FontSize',12);
title('Cold flow varying','FontSize',12);
legend('Location','best','FontSize',10); grid on; box on; xlim([0 2000]);
set(gca,'FontSize',11);

subplot(1,2,2); hold on;
plot(Q_sweep, Th_cc_hot.psi, '-',  'Color',blue,  'LineWidth',lw_th, ...
     'DisplayName','Counter (theory)');
plot(Q_sweep, Th_ph_hot.psi, '--', 'Color',coral, 'LineWidth',lw_th, ...
     'DisplayName','Parallel (theory)');
scatter(x_hot_C, E_CH.psi, 80, blue,  'o','filled','DisplayName','Counter (exp.)');
scatter(x_hot_P, E_PH.psi, 80, coral, 's','filled','DisplayName','Parallel (exp.)');
xlabel('Hot flow rate (l/h)','FontSize',12);
ylabel('\psi (%)','FontSize',12);
title('Hot flow varying','FontSize',12);
legend('Location','best','FontSize',10); grid on; box on; xlim([0 2000]);
set(gca,'FontSize',11);

%% =========================================================================
%  SECTION 10 — PRINT SUMMARY TABLE
% =========================================================================
fprintf('\n%s\n', repmat('=',1,125));
fprintf('  SUMMARY TABLE — experimental vs theoretical (at matched operating points)\n');
fprintf('%s\n', repmat('=',1,125));
fprintf('%-22s  %-7s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s\n', ...
    'Dataset','Q(l/h)','Qh_e','Qh_t','I_e','I_t','Sg_e','Sg_t','eps_e','eps_t','psi_e','psi_t');
fprintf('%s\n', repmat('-',1,125));

% Compute theory at exact points for table
all_data   = {data_CC, data_PC, data_CH, data_PH};
all_exp    = {E_CC, E_PC, E_CH, E_PH};
all_config = {'counter','parallel','counter','parallel'};
all_sweep  = {'cold','cold','hot','hot'};
all_labels = {'Counter cold-var','Parallel cold-var','Counter hot-var','Parallel hot-var'};

for d = 1:4
    data_d  = all_data{d};
    E_d     = all_exp{d};
    cfg     = all_config{d};
    sw      = all_sweep{d};
    lbl     = all_labels{d};
    for i = 1:size(data_d,1)
        r = compute_theory_point(data_d(i,1),data_d(i,2),data_d(i,3),data_d(i,5),cfg,T0,...
            rho_w,cp_w,mu_w,lam_w,Pr_w,di,de,nt,L,k_w,Sm,A_sh,Dh,Si,Se_a);
        if strcmp(sw,'cold')
            q_print = data_d(i,2);
        else
            q_print = data_d(i,1);
        end
        fprintf('%-22s  %-7.0f  %-9.1f  %-9.1f  %-9.3f  %-9.3f  %-9.1f  %-9.1f  %-9.2f  %-9.2f  %-9.2f  %-9.2f\n', ...
            lbl, q_print, ...
            E_d.Qh(i),  r.Qh, ...
            E_d.Idest(i), r.Idest, ...
            E_d.Sgen(i),  r.Sgen, ...
            E_d.eps(i),   r.eps, ...
            E_d.psi(i),   r.psi);
        lbl = '';
    end
    fprintf('%s\n', repmat('-',1,125));
end

fprintf('\nNotes:\n');
fprintf('  Qh  = heat released by hot fluid [W]\n');
fprintf('  I   = exergy destruction = dXh - dXc [W]\n');
fprintf('  Sg  = entropy generation × 1000 [mW/K]\n');
fprintf('  eps = 1st-law effectiveness [%%]\n');
fprintf('  psi = I/dXh × 100 [%%]  (fraction of hot-side exergy destroyed)\n');
fprintf('        Cold inlet temps 6.4–8°C < T0=10°C, so dXc<0 in all tests.\n');
fprintf('        Classical psi=dXc/dXh is undefined here; I/dXh is used instead.\n');
fprintf('        Lower psi = less exergy wasted → counter-flow is more efficient.\n');
fprintf('\nTheoretical curves span 10–2000 l/h using mean inlet conditions of each dataset.\n');
fprintf('All done. Figures 1–4 generated.\n');