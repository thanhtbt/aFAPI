% Effect of parameters alpha and p on performance of aFAPI
clear;clc; close all
addpath(genpath('subspace_trackers\'))

%% Parameters
n = 50;
r = 5;
T = 1001;
time_varying_factor = 1e-3*ones(1,T);
beta    = 0.99; % forgetting factor
alpha   = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99] ;  % alpha-divergence
lp      = [1 1.25 1.5 1.75 2];

% Generating data
% true data
disp('+ Data Generating ...')    
[X,U_tr] = data_generator(n,T,r,time_varying_factor);
% noise
fac_noise = 1;
epsilon   = 0.2;  % for non-gaussian noise
sigma_n   = 10;
mu_n      = 10;

Noise     = (1-epsilon)*randn(n,T) + epsilon*(randn(n,T));
outlier_1 = (1-epsilon)*randn(n,1) + epsilon*0.5*sigma_n*(randn(n,1) + mu_n);  
outlier_2 = (1-epsilon)*randn(n,1) + epsilon*sigma_n*(randn(n,1) + mu_n);  
outlier_3 = (1-epsilon)*randn(n,1) + epsilon*2*sigma_n*(randn(n,1) + mu_n);  
Noise(:,400)  = outlier_2;
Noise(:,600)  = outlier_3;
Noise(:,800)  = outlier_3;


X_noise   = X + Noise;

disp('+ Processing ...')
idx = 1;
for ii = 1 : length(alpha)
    for jj = 1 : length(lp)
        fprintf('+ Case %d: (alpha = %0.1f | lp = %0.1f) \n',idx,alpha(ii),lp(jj))
        [~, eta_RFAPI_ii,rho_RFAPI_ii]  = alpha_FAPI(X_noise,beta,alpha(ii),U_tr,lp(jj));
        eta_RFAPI{ii,jj}  = eta_RFAPI_ii;
        rho_RFAPI{ii,jj}  = rho_RFAPI_ii;
        idx = idx + 1;
    end
end

%% PLOT
disp('Plotting ....')

fig1 = figure; hold on;
for ii = 1:7
    semilogy(eta_RFAPI{ii,5},'LineWidth',1.5);
end
semilogy(eta_RFAPI{8,5},'-b','LineWidth',2);
semilogy(eta_RFAPI{9,5},'-g','LineWidth',2);
semilogy(eta_RFAPI{10,5},'-r','LineWidth',2);
axis([0 1000 1e-3 10 ]);

leg = legend('$\alpha=0.1$','$\alpha=0.2$','$\alpha=0.3$','$\alpha=0.4$','$\alpha=0.5$', ...
             '$\alpha=0.6$','$\alpha=0.7$','$\alpha=0.8$','$\alpha=0.9$','$\alpha=0.99$')
set(leg,'interpreter','latex');
title('p = 2 (standard)')
ylabel('SEP','interpreter','latex','FontSize',13,'FontName','Times New Roman'); 
xlabel('Data Samples','interpreter','latex','FontSize',13,'FontName','Times New Roman');
h=gca;
set(gca, 'YScale', 'log','FontSize',24)
set(h,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
set(h,'Xtick',0:200:T,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
grid on;
box on;
set(fig1, 'units', 'inches', 'position', [0.5 0.5 9 7]);


fig2 = figure; hold on;
for ii = 1:5
    semilogy(eta_RFAPI{7,ii},'LineWidth',1.5);
end

axis([0 1000 1e-3 10 ]);
ylabel('SEP','interpreter','latex','FontSize',13,'FontName','Times New Roman'); 
xlabel('Data Samples','interpreter','latex','FontSize',13,'FontName','Times New Roman');
h=gca;
set(gca, 'YScale', 'log','FontSize',24)
set(h,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
set(h,'Xtick',0:200:T,'FontSize',24,'XGrid','on','YGrid','on','FontName','Times New Roman');
grid on;
box on;
set(fig2, 'units', 'inches', 'position', [0.5 0.5 9 7]);
leg2 = legend('$p=1$','$p=1.25$','$p=1.5$','$p=1.75$','$p=2$');
set(leg2,'interpreter','latex');
title('$\alpha = 0.7$','interpreter','latex')
 