% Adapted by Frank DeStefano
% Adapted by William Ojemann
% BioE 485 Final Project
% 5/24/2020
% This file is used to run the minimize the objective function and fit the
% parameters to data from Snooke et al.
% Based on code from Dr. Thomas
 close all; clc;
% this is original high density lifetime data
data = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        18.0890E-12    1.8626
        29.8906E-12    1.3289
        39.2560E-12    1.3382];
  % this is standard deviation high density data
%     data = [5.981995393982141E-12   1.117118131404848
%   11.879121014737230E-12   1.258281889845138
%   18.041235762195786E-12   1.359373827229811
%   29.853021931033279E-12   1.360952954674168
%   39.127711940341570E-12   1.358182170871337];

% Initial guesses


kc = 10;
% ks = .25;
k01 = 1e-12; % rate constant from sun et al
D = .01e-8;
kappa = [0.00131211798645895];
eta = [1.51813939524036e-10];
ks = [0.255082466895303];
%kc = [1.69873006371107];

% Parameter Fitting
% choose either new parameters to fit or go with last estimates
%guesses = [kappa, eta, ks, kc, k01, D];
 %guesses = [estimates];
bestsofar = [0.00124402839199090,2.55295462419621e-10,0.111384250721051,1.56374425143498,1,8.03092771313965e-11];
guesses = bestsofar;
figure
[estimates,J] = fminsearch(@objj,guesses);
disp(estimates)
tmean = variedIClifetime(estimates,.58,data(:,1));




% good est = [0.00283356626470898,3.73283447342488e-10,0.0645436399125567,2.58974669653954,6.49143167265349e-14,1.31852729200664e-09]
bestsofar = [0.00124402839199090,2.55295462419621e-10,0.111384250721051,1.56374425143498,1.18318602769146e-12,8.03092771313965e-11];
 