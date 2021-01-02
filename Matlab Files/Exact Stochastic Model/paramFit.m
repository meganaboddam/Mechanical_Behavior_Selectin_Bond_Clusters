% Jasmeet Khera (David, Frank, Megana, Will)
% Adapted by William Ojemann
% BioE 485 Final Project
% 5/24/2020
% This file is used to run the minimize the objective function and fit the
% parameters to data from Snooke et al.
% Based on code from Dr. Thomas
clear all; close all; clc;
global data
data = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        18.0890E-12    1.8626
        29.8906E-12    1.3289
        39.2560E-12    1.3382];

% Initial guesses
N0 = 10;
kappa = 10;
k0b = 1;
eta = 1E-8;
kc = .001;
ks = 300;

% Parameter Fitting
guesses = [kappa, k0b,eta,kc,ks,N0];
figure
[estimates,J] = fminsearch(@objective,guesses);
disp(estimates)
[forces,tmean] = lifetimesExact(estimates);
y2 = interp1(forces,tmean,data(:,1));
figure
plot(data(:,1),y2,data(:,1),data(:,2),'o');

% Objective function for evaluating parameter fitting
function J = objective(guesses)
disp(guesses)
global data
[forces, tmean] = lifetimesExact(guesses);
y2 = interp1(forces,tmean,data(:,1));
plot(data(:,1),y2,'--',data(:,1),data(:,2),'o'); drawnow
J = sum((data(:,2)-y2).^2);
end

