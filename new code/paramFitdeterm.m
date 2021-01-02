% Adapted by Megana Boddam
% Adapted by Frank DeStefano
% Adapted by William Ojemann
% BioE 485 Final Project
% 5/24/2020
% This file is used to create the force sensitive and insensitive
% verification figures
% This file is used to run the minimize the objective function and fit the
% parameters to data from Snooke et al.
% Based on code from Dr. Thomas
 close all; clear all;

format long;
 
%                           Data format - [ force, lifetime]
%
%                         --High Density (Cluster) data -- 
% Mean Lifetime
% HDlife = [5.9519E-12  1.3129
%         12.0350E-12    1.3329
%         18.0890E-12    1.8626
%         29.8906E-12    1.3289
%         39.2560E-12    1.3382];
% 
% 
% % Standard Deviation
% HDstd= [5.981995393982141E-12   1.117118131404848
%            11.879121014737230E-12   1.258281889845138
%            18.041235762195786E-12   1.359373827229811
%            29.853021931033279E-12   1.360952954674168
%            39.127711940341570E-12   1.358182170871337];
       
% Data after removing the outlier (3rd point)
HDlife_noOutlier = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        29.8906E-12    1.3289
        39.2560E-12    1.3382];
HDstd_noOutlier = [5.981995393982141E-12   1.117118131404848
           11.879121014737230E-12   1.258281889845138
           29.853021931033279E-12   1.360952954674168
           39.127711940341570E-12   1.358182170871337];       
    
    
%                           --Low Density (single) data --
%Mean Lifetime
LDlife = [6.04805000000000e-12,0.433328000000000;
         1.20022000000000e-11,0.806065000000000;
         1.80241000000000e-11,0.912796000000000;
         2.98636000000000e-11,0.799691000000000;
         3.91329000000000e-11,0.296444000000000];

%Standard Deviation
LDstd= [6.09917e-12	0.635847; 
         12.2032e-12 1.06711;
         18.1886e-12	1.16884;
         30.1231e-12	1.01431;
         39.2188e-12	0.435177];

% Data after removing the outlier (3rd point)
LDlife_noOutlier = [6.04805000000000e-12,0.433328000000000;
         1.20022000000000e-11,0.806065000000000;
         2.98636000000000e-11,0.799691000000000;
         3.91329000000000e-11,0.296444000000000];

LDstd_noOutlier = [6.09917e-12	0.635847; 
         12.2032e-12 1.06711;
         30.1231e-12	1.01431;
         39.2188e-12	0.435177];

% figure 1
% use initial guesses to find optimal parameters for low density cluster to show
% catch behavior. 

% Initial guesses from general guesses after reading various articles
D = 1.319e-9; 
eta = 3.73e-10;
kc = 2.59;
ks = 6.45e-2;
k01 = 6.491e-14;
kappa = 2.833e-3;
guesses1 = [kappa, eta, ks, kc, k01, D];

f = linspace(0, 40e-12, 50);
[time1,~] = variedIClifetime2(guesses1, 0.25, f);

figure(1);
plot(f, time1);
ylim([0 2]);
xlim([0 40e-12]);
title("Low Density Catch Bond Cluster Lifetimes Over Force", 'FontSize', 12);
xlabel("Force (pN)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);

%% figure 2
% set eta and characteristic length to 0 to see if results end up being
% force insensitive. Rest of initial conditions same as in figure 1.
format long;
D = 0; eta = 0;
kc = 2.59;
ks = 6.45e-2;
k01 = 6.491e-14;
kappa = 2.833e-3;
guesses2 = [kappa, eta, ks, kc, k01, D];

f = linspace(0, 40e-12, 50);
[time2,~] = variedIClifetime2(guesses2, 0.25, f);

figure(2);
plot(f, time2, 'LineWidth', 2);
ylim([0 2]);
title("Low Density Catch Bond Cluster Lifetimes Over Force", 'FontSize', 12);
xlabel("Force (pN)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);

%% figures 3-6
% Validation. Parameter Fitting of Model to High Density Data
format long;

kc = 1;
ks = .25;
k01 = 4;
D = .01e-6;
eta = 1.51813939524036e-10;
kappa = 0.00131211798645895;
guesses3 = [kappa, eta, ks, kc, k01, D];

[estimates_noOutlier_HDfit, J] = fminsearch(@objjNoOutlierHDFit, guesses3);
disp(estimates_noOutlier_HDfit);
disp(J);

[time,~] = variedIClifetime2(estimates_noOutlier_HDfit, 0.58, HDlife_noOutlier(:,1));
[~,std] = variedIClifetime2(estimates_noOutlier_HDfit, 0.58, HDstd_noOutlier(:,1)); 

figure(3);
plot(HDlife_noOutlier(:,1), time, '*', 'LineWidth', 2); hold on;
plot(HDlife_noOutlier(:,1), HDlife_noOutlier(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("Model Fit to High Density Mean Data (excluding outlier)", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "mean HD data", 'FontSize', 12);

figure(4);
plot(HDstd_noOutlier(:,1), std,'*', 'LineWidth', 2); hold on;  
plot(HDstd_noOutlier(:,1),HDstd_noOutlier(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("High Density Standard Deviation Data (excluding outlier)", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "stdev HD data", 'FontSize', 12);

[time1,~] = variedIClifetime2(estimates_noOutlier_HDfit, 0.25, LDlife(:,1));
[~,std2] = variedIClifetime2(estimates_noOutlier_HDfit, 0.25, LDstd(:,1)); 

figure(5);
plot(LDlife(:,1), time1, '*', 'LineWidth', 2); hold on;
plot(LDlife(:,1), LDlife(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("High Density Model Parameters Fit to Low Density Data (excluding outlier)", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "mean LD data", 'FontSize', 12);

figure(6);
plot(LDstd(:,1), std2,'*', 'LineWidth', 2); hold on;  
plot(LDstd(:,1),LDstd(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("Low Density Standard Deviation (excluding outlier)", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "stdev LD data", 'FontSize', 12);

%% figures 7-10
% Validation. Parameter Fitting of Model to Low Density Data
kc = 1;
ks = .25;
k01 = 4;
D = .01e-6;
eta = 1.51813939524036e-10;
kappa = 0.00131211798645895;
guesses4 = [kappa, eta, ks, kc, k01, D];

[estimates_LDfit, J2] = fminsearch(@objjLDFit, guesses4);
disp(estimates_LDfit);
disp(J2);

[time3,~] = variedIClifetime2(estimates_LDfit, 0.25, LDlife(:,1));
[~,std3] = variedIClifetime2(estimates_LDfit, 0.25, LDstd(:,1)); 

figure(7);
plot(LDlife(:,1), time3, '*', 'LineWidth', 2); hold on;
plot(LDlife(:,1), LDlife(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("Model Fit to Low Density Data", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "mean LD data", 'FontSize', 12);

figure(8);
plot(LDstd(:,1), std3,'*', 'LineWidth', 2); hold on;  
plot(LDstd(:,1),LDstd(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("Low Density Standard Deviation Data", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "stdev LD data", 'FontSize', 12);

[time4,~] = variedIClifetime2(estimates_LDfit, 0.58, HDlife_noOutlier(:,1));
[~,std4] = variedIClifetime2(estimates_LDfit, 0.58, HDstd_noOutlier(:,1)); 

figure(9);
plot(HDlife_noOutlier(:,1), time4, '*', 'LineWidth', 2); hold on;
plot(HDlife_noOutlier(:,1), HDlife_noOutlier(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("Low Density Model Parameters Fit to High Density Data (excluding outlier)", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "mean HD data", 'FontSize', 12);

figure(10);
plot(HDstd_noOutlier(:,1), std4,'*', 'LineWidth', 2); hold on;  
plot(HDstd_noOutlier(:,1),HDstd_noOutlier(:,2),'o', 'LineWidth', 2); hold off;
ylim([0 2]);
title("High Density Standard Deviation Data (excluding outlier)", 'FontSize', 12);
xlabel("Force (N)", 'FontSize', 12);
ylabel("Bond Lifetime (s)", 'FontSize', 12);
legend("model", "stdev HD data", 'FontSize', 12);


