% Demo Subspace Tracking
clear;clc; close all
addpath(genpath('subspace_trackers\'))

%% Experimental parameters
n_exp = 10;       % Number of independent runs
n     = 50;       % Data dimension
r     = 5;        % Target rank
T     = 1000;     % Data samples
beta  = 0.99;     % forgetting factor
alpha = 0.9;      % alpha-divergence
time_varying_factor = 1e-3*ones(1,T); % time-varying


%% Performance metrics
eta_OPAST  = zeros(1,T); rho_OPAST  = zeros(1,T);
eta_aFAPI  = zeros(1,T); rho_aFAPI  = zeros(1,T);
eta_FAPI   = zeros(1,T); rho_FAPI   = zeros(1,T);
eta_RPAST  = zeros(1,T); rho_RPAST  = zeros(1,T);
eta_YAST   = zeros(1,T); rho_YAST   = zeros(1,T);
eta_LORAF  = zeros(1,T); rho_LORAF  = zeros(1,T);
eta_TRPAST = zeros(1,T); rho_TRPAST = zeros(1,T);
eta_aFAPI  = zeros(1,T); rho_aFAPI  = zeros(1,T);


for ii = 1 : n_exp
    disp('-------------------------')
    fprintf('Run %d/%d \n',ii,n_exp)
    
    disp('+ Data Generating ...')
    % Generating data
    % true data
    [X,U_tr]  = data_generator(n,T,r,time_varying_factor);
    % Contaminated Mixture Noises
    fac_noise = 1;
    epsilon   = 0.2;
    sigma_n   = 10;
    mu_n      = 10;
    Noise     = (1-epsilon)*randn(n,T) + epsilon*(randn(n,T));
    
    % Abrupt changes
    outlier_1 = (1-epsilon)*randn(n,1) + epsilon*0.5*sigma_n*(randn(n,1) + mu_n);
    outlier_2 = (1-epsilon)*randn(n,1) + epsilon*sigma_n*(randn(n,1) + mu_n);
    outlier_3 = (1-epsilon)*randn(n,1) + epsilon*2*sigma_n*(randn(n,1) + mu_n);
    Noise(:,400)  = outlier_2; % 1.5*sigma_n*rand(n,1);
    Noise(:,600)  = outlier_3; % 1.5*sigma_n*rand(n,1);
    Noise(:,800)  = outlier_3; % 1.5*sigma_n*rand(n,1);
    
    X_noise   = X + Noise;
    
    %% Main Program
    disp('+ Processing ...')
    disp('    - aFAPI')
    [~, eta_aFAPI_ii,rho_aFAPI_ii]     = alpha_FAPI(X_noise,beta,alpha,U_tr,1.5);
    disp('    - FAPI')
    [~,eta_FAPI_ii,rho_FAPI_ii]        = FAPI(X_noise,beta,U_tr);
    disp('    - OPAST')
    [~,eta_OPAST_ii,rho_OPAST_ii]      = OPAST(X_noise,beta,U_tr);
    disp('    - YAST')
    [~,eta_YAST_ii,rho_YAST_ii]        = GYAST(X_noise,beta,U_tr);
    disp('    - RPAST')
    [~,eta_RPAST_ii,rho_RPAST_ii]      = RPAST(X_noise,beta,U_tr);
    disp('    - RRPAST')
    [~, eta_TRPAST_ii,rho_TRPAST_ii]   = TRPAST(X_noise,beta,alpha,U_tr);
    disp('    - LORAF')
    [~, eta_LORAF_ii,rho_LORAF_ii]     = LORAF(X_noise,beta,U_tr);
    
    eta_aFAPI   = eta_aFAPI  + eta_aFAPI_ii;
    eta_FAPI    = eta_FAPI   + eta_FAPI_ii;
    eta_OPAST   = eta_OPAST  + eta_OPAST_ii;
    eta_RPAST   = eta_RPAST  + eta_RPAST_ii;
    eta_YAST    = eta_YAST   + eta_YAST_ii;
    eta_TRPAST  = eta_TRPAST + eta_TRPAST_ii;
    eta_LORAF   = eta_LORAF  + eta_LORAF_ii;
    
end
eta_aFAPI   = eta_aFAPI/n_exp;
eta_FAPI    = eta_FAPI/n_exp;
eta_RPAST   = eta_RPAST/n_exp;
eta_OPAST   = eta_OPAST/n_exp;
eta_YAST    = eta_YAST/n_exp;
eta_TRPAST  = eta_TRPAST/n_exp;
eta_LORAF   = eta_LORAF /n_exp;

%% PLOT
disp('+ Plotting ....')

makerSize   = 11;
numbMarkers = 500;
LineWidth   = 2;

color      = get(groot,'DefaultAxesColorOrder');
red_o      = [1,0,0];
blue_o     = [0, 0, 1];
magenta_0  = [1 0 1];
gree_o     = [0, 0.5, 0];
black_o    = [0.25, 0.25, 0.25];
blue_n     = color(1,:);
oran_n     = color(2,:);
yell_n     = color(3,:);
viol_n     = color(4,:);
gree_n     = color(5,:);
lblu_n     = color(6,:);
brow_n     = color(7,:);
lbrow_n    = [0.5350    0.580    0.2840];

k = 1; TS = 100;
fig = figure; hold on;

t1  = semilogy(1:k:T,eta_FAPI(1:k:T),'--b','Color',viol_n,'LineWidth',LineWidth);
t11 = semilogy(1:TS:T,eta_FAPI(1:TS:T),'marker','x','markersize',makerSize,...
    'linestyle','none','color',viol_n,'LineWidth',LineWidth);
t12 =  semilogy(1,eta_FAPI(1),'marker','x','markersize',makerSize,...
    'linestyle','-.','color',viol_n,'LineWidth',LineWidth);

t2  = semilogy(1:k:T,eta_YAST(1:k:T),'-','color',gree_o,'LineWidth',LineWidth);
t21 = semilogy(1:TS:T,eta_YAST(1:TS:T),'marker','+','markersize',makerSize,...
    'linestyle','none','color',gree_o,'LineWidth',LineWidth);
t22 =  semilogy(1,eta_YAST(1),'marker','+','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

t3  = semilogy(1:k:T,eta_LORAF(1:k:T),'--','color',magenta_0,'LineWidth',1.5);
t31 = semilogy(1:TS:T,eta_LORAF(1:TS:T),'marker','o','markersize',makerSize,...
    'linestyle','none','color',magenta_0,'LineWidth',LineWidth);
t32 =  semilogy(1,eta_LORAF(1),'marker','o','markersize',makerSize,...
    'linestyle','-.','color',magenta_0,'LineWidth',LineWidth);

t4 = semilogy(1:k:T,eta_OPAST(1:k:T),'-','color',blue_o,'LineWidth',1.5);
t41 = semilogy(1:TS:T,eta_OPAST(1:TS:T),'marker','^','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
t42 =  semilogy(1,eta_OPAST(1),'marker','^','markersize',makerSize,...
    'linestyle','-.','color',blue_o,'LineWidth',LineWidth);

t6  = semilogy(1:k:T,eta_TRPAST(1:k:T),'-m','color',black_o,'LineWidth',LineWidth);
t61 =  semilogy(1:TS:T,eta_TRPAST(1:TS:T),'marker','d','markersize',makerSize,...
    'linestyle','none','color',black_o,'LineWidth',LineWidth);
t62 =  semilogy(1,eta_TRPAST(1),'marker','d','markersize',makerSize,...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);

t7  = semilogy(1:k:T,eta_aFAPI(1:k:T),'-','color',red_o,'LineWidth',2.5);
t71 =  semilogy(1:TS:T,eta_aFAPI(1:TS:T),'marker','p','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
t72 =  semilogy(1,eta_aFAPI(1),'marker','p','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

hold off

ylabel('SEP','interpreter','latex','FontSize',13,'FontName','Times New Roman');
xlabel('Data Samples','interpreter','latex','FontSize',13,'FontName','Times New Roman');
leg = legend([t12 t32 t22 t42  t62 t72],...
    '\texttt{FAPI}','\texttt{LORAF}','\texttt{YAST}','\texttt{OPAST}','\texttt{TRPAST}','$\alpha$\texttt{FAPI(Proposed)}');
set(leg,'Location','NorthEast','Interpreter','latex',...
    'FontSize',22,'NumColumns',2,'EdgeColor', 'none');

axis([0 T 1e-3 10]);
h=gca;
set(gca, 'YScale', 'log','FontSize',24)
set(h,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
set(h,'Xtick',0:round(T/5):T,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 9 7]);

