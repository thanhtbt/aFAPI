% Uniform Linear Array(ULA) Model
% DOA tracking
clear;close all; clc
addpath(genpath('subspace_trackers\'))

%% Experimental parameters
n      = 20;        % number of sensors
p      = 3;         % number of sources
T      = 1000;      % number of data samples
beta   = 0.99;      % forgetting factor

%% Data Generating
disp('Data generating ...')
% create original signals
S = (randn(p,T) + 1i*randn(p,T))/sqrt(2);
% construct matrix A
theta1 = (20:20/(T-1):40)'*pi/180;
theta2 = (10:-10/(T-1):0)'*pi/180;
theta3 = -(10*pi/180)*ones(T,1);
theta  = [theta1 theta2 theta3];

lamb = 1; d = lamb/2;

X = zeros(n,T);
for tt = 1 : T
    w  = 2*pi*(d/lamb)*sin(theta(tt,:));
    A  = exp(1i*(0:n-1)'*w);     % A changes slowly
    X(:,tt) = A*S(:,tt);
end
Omega = 0.5*sin(theta);

%% Noise
% Background noises
snr     = 0; % dB
std_brt = 10^(-snr/20);   
G       = std_brt * (randn(n, T)+ 1i * randn(n, T))/sqrt(2); 
X       = X + G;
% Abrupt noises
epsilon = 0.2;
t       = 500:599;
Nt      = (1-epsilon)*randn(n,length(t)) + epsilon*5*(randn(n,length(t))+10);
X(:,t)  = X(:,t)+ Nt;

%% DOA Tracking
disp('Tracking ...')

disp(' - FAPI')
estOmg_FAPI = DOA_FAPI(X,beta,p);

disp(' - OPAST')
estOmg_OPAST = DOA_OPAST(X,beta,p);

disp(' - LORAF')
estOmg_LORAF = DOA_LORAF(X,beta,p);

disp(' - TRPAST')
alpha  = 0.9;
estOmg_TRPAST = DOA_TRPAST(X,beta,alpha,p);

disp(' - alpha FAPI')
alpha  = 0.9;
estOmg_alpha_FAPI = DOA_alpha_FAPI(X,beta,p,alpha);


%% PLOT

disp('Plotting ....')

makerSize   = 11;
numbMarkers = 500;
LineWidth   = 2;

color      = get(groot,'DefaultAxesColorOrder');
red_o      = [1,0,0];
blue_o     = [0, 0, 1];
magenta_0  = [1 0 1];
gree_o     = [0, 0.5, 0];
black_o    = [0.25, 0.25, 0.25];
lbrow_n    = [0.5350 0.580 0.2840];
blue_n     = color(1,:);
oran_n     = color(2,:);
yell_n     = color(3,:);
viol_n     = color(4,:);
gree_n     = color(5,:);
lblu_n     = color(6,:);
brow_n     = color(7,:);

k = 1; TS = 100;

fig = figure; hold on;

% Source 1
t0 = plot(1:k:T,Omega(1:k:T,1),'-','Color','k','LineWidth',3); 
t01 = plot(1:TS:T,Omega(1:TS:T,1),'markersize',makerSize,...
   'linestyle','none','color','k','LineWidth',3);
t02 =  plot(1,Omega(1,1),'markersize',makerSize,...
   'linestyle','-','color','k','LineWidth',3);


t1 = plot(1:k:T,estOmg_FAPI(3,1:k:T),'-','Color',viol_n,'LineWidth',LineWidth); 
t11 = plot(1:TS:T,estOmg_FAPI(3,1:TS:T),'marker','x','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',LineWidth);
t12 =  plot(1,estOmg_FAPI(3,1),'marker','x','markersize',makerSize,...
   'linestyle','-','color',viol_n,'LineWidth',LineWidth);


t2 = plot(1:k:T,estOmg_OPAST(3,1:k:T),'-','Color','g','LineWidth',LineWidth); 
t21 = plot(1:TS:T,estOmg_OPAST(3,1:TS:T),'marker','s','markersize',makerSize,...
   'linestyle','none','color','g','LineWidth',LineWidth);
t22 =  plot(1,estOmg_OPAST(3,1),'marker','s','markersize',makerSize,...
   'linestyle','-','color','g','LineWidth',LineWidth);


t3 = plot(1:k:T,estOmg_LORAF(3,1:k:T),'-','Color','b','LineWidth',LineWidth); 
t31 = plot(1:TS:T,estOmg_LORAF(3,1:TS:T),'marker','^','markersize',makerSize,...
   'linestyle','none','color','b','LineWidth',LineWidth);
t32 =  plot(1,estOmg_LORAF(3,1),'marker','^','markersize',makerSize,...
   'linestyle','-','color','b','LineWidth',LineWidth);


t4 = plot(1:k:T,estOmg_TRPAST(3,1:k:T),'-','Color',gree_o,'LineWidth',LineWidth); 
t41 = plot(1:TS:T,estOmg_TRPAST(3,1:TS:T),'marker','d','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
t42 =  plot(1,estOmg_TRPAST(3,1),'marker','d','markersize',makerSize,...
   'linestyle','-','color',gree_o,'LineWidth',LineWidth);

t5 = plot(1:k:T,estOmg_alpha_FAPI(3,1:k:T),'-','Color','r','LineWidth',LineWidth); 
t51 = plot(1:TS:T,estOmg_alpha_FAPI(3,1:TS:T),'marker','p','markersize',makerSize,...
   'linestyle','none','color','r','LineWidth',LineWidth);
t52 =  plot(1,estOmg_alpha_FAPI(3,1),'marker','p','markersize',makerSize,...
   'linestyle','-','color','r','LineWidth',LineWidth);

% Source 2
t6 = plot(1:k:T,Omega(1:k:T,2),'-','Color','k','LineWidth',3); 
t61 = plot(1:TS:T,Omega(1:TS:T,2),'markersize',makerSize,...
   'linestyle','none','color',brow_n,'LineWidth',3);
t62 =  plot(1,Omega(1,1),'markersize',makerSize,...
   'linestyle','-','color',brow_n,'LineWidth',3);

t7 = plot(1:k:T,estOmg_FAPI(2,1:k:T),'-','Color',viol_n,'LineWidth',LineWidth); 
t71 = plot(1:TS:T,estOmg_FAPI(2,1:TS:T),'marker','x','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',LineWidth);
t72 =  plot(1,estOmg_FAPI(2,1),'marker','x','markersize',makerSize,...
   'linestyle','-','color',viol_n,'LineWidth',LineWidth);

t8 = plot(1:k:T,estOmg_OPAST(2,1:k:T),'-','Color','g','LineWidth',LineWidth); 
t81 = plot(1:TS:T,estOmg_OPAST(2,1:TS:T),'marker','s','markersize',makerSize,...
   'linestyle','none','color','g','LineWidth',LineWidth);
t82 =  plot(1,estOmg_OPAST(2,1),'marker','s','markersize',makerSize,...
   'linestyle','-','color','g','LineWidth',LineWidth);

t9 = plot(1:k:T,estOmg_LORAF(2,1:k:T),'-','Color','b','LineWidth',LineWidth); 
t91 = plot(1:TS:T,estOmg_LORAF(2,1:TS:T),'marker','^','markersize',makerSize,...
   'linestyle','none','color','b','LineWidth',LineWidth);
t92 =  plot(1,estOmg_LORAF(2,1),'marker','^','markersize',makerSize,...
   'linestyle','-','color','b','LineWidth',LineWidth);

t10_ = plot(1:k:T,estOmg_TRPAST(2,1:k:T),'-','Color',gree_o,'LineWidth',LineWidth); 
t101 = plot(1:TS:T,estOmg_TRPAST(2,1:TS:T),'marker','d','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
t102 =  plot(1,estOmg_TRPAST(2,1),'marker','d','markersize',makerSize,...
   'linestyle','-','color',gree_o,'LineWidth',LineWidth);

t11_ = plot(1:k:T,estOmg_alpha_FAPI(2,1:k:T),'-','Color','r','LineWidth',LineWidth); 
t111 = plot(1:TS:T,estOmg_alpha_FAPI(2,1:TS:T),'marker','p','markersize',makerSize,...
   'linestyle','none','color','r','LineWidth',LineWidth);
t112 =  plot(1,estOmg_alpha_FAPI(2,1),'marker','p','markersize',makerSize,...
   'linestyle','-','color','r','LineWidth',LineWidth);

% Source 3
t120 = plot(1:k:T,Omega(1:k:T,3),'-','Color','k','LineWidth',3); 
t121 = plot(1:TS:T,Omega(1:TS:T,3),'markersize',makerSize,...
   'linestyle','none','color',brow_n,'LineWidth',3);
t122 =  plot(1,Omega(1,1),'markersize',makerSize,...
   'linestyle','--','color',brow_n,'LineWidth',3);

t130 = plot(1:k:T,estOmg_FAPI(1,1:k:T),'-','Color',viol_n,'LineWidth',LineWidth); 
t131 = plot(1:TS:T,estOmg_FAPI(1,1:TS:T),'marker','x','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',LineWidth);
t132 =  plot(1,estOmg_FAPI(1,1),'marker','x','markersize',makerSize,...
   'linestyle','-','color',viol_n,'LineWidth',LineWidth);

t140 = plot(1:k:T,estOmg_OPAST(1,1:k:T),'-','Color','g','LineWidth',LineWidth); 
t141 = plot(1:TS:T,estOmg_OPAST(1,1:TS:T),'marker','s','markersize',makerSize,...
   'linestyle','none','color','g','LineWidth',LineWidth);
t142 =  plot(1,estOmg_OPAST(1,1),'marker','s','markersize',makerSize,...
   'linestyle','-','color','g','LineWidth',LineWidth);

t150 = plot(1:k:T,estOmg_LORAF(1,1:k:T),'-','Color','b','LineWidth',LineWidth); 
t151 = plot(1:TS:T,estOmg_LORAF(1,1:TS:T),'marker','^','markersize',makerSize,...
   'linestyle','none','color','b','LineWidth',LineWidth);
t152 =  plot(1,estOmg_LORAF(1,1),'marker','^','markersize',makerSize,...
   'linestyle','-','color','b','LineWidth',LineWidth);

t160 = plot(1:k:T,estOmg_TRPAST(1,1:k:T),'-','Color',gree_o,'LineWidth',LineWidth); 
t161 = plot(1:TS:T,estOmg_TRPAST(1,1:TS:T),'marker','d','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
t162 =  plot(1,estOmg_TRPAST(2,1),'marker','d','markersize',makerSize,...
   'linestyle','-','color',gree_o,'LineWidth',LineWidth);

t170 = plot(1:k:T,estOmg_alpha_FAPI(1,1:k:T),'-','Color','r','LineWidth',LineWidth); 
t171 = plot(1:TS:T,estOmg_alpha_FAPI(1,1:TS:T),'marker','p','markersize',makerSize,...
   'linestyle','none','color','r','LineWidth',LineWidth);
t172 =  plot(1,estOmg_alpha_FAPI(1,1),'marker','p','markersize',makerSize,...
   'linestyle','-','color','r','LineWidth',LineWidth);

ylabel('Angle Frequency','interpreter','latex','FontName','Times New Roman'); 
xlabel('Data Sample','interpreter','latex','FontName','Times New Roman');
leg = legend([ t12 t22 t32 t42 t52 t02],...
    '\texttt{FAPI}','\texttt{OPAST}','\texttt{LORAF}',...
    '\texttt{TRPAST}','$\alpha$\texttt{FAPI(Proposed)}','\texttt{Ground truth}');
set(leg,'Location','NorthWest','Interpreter','latex',...
    'FontSize',18,'NumColumns',2,'EdgeColor', 'none');
title(sprintf('DOA tracking'),'interpreter','latex','FontSize',14);

h = gca;
set(gca,'FontSize',24)
set(h,'Xtick',0:200:1000,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
grid on;
box on;
axis([0 1000 -0.2 0.4])
set(fig, 'units', 'inches', 'position', [0.5 0.5 10 7]);




