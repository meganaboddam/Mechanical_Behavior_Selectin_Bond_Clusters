% William Ojemann
% BioE 485 Final Project
% General Script
clear all; close all; clc;
% Parameters are chosen based off of Sun et al.
tstop = 200;
N0 = 5;
kappa = 10;
% these parameters are estimates from Novikova et al.
k0b = 1;
eta = 1E-8;
kc = .001;
ks = 300;
sample = 100;
forces = linspace(1,50)*1E-12;
tmean = zeros(1,length(forces));
for i = 1:length(forces)
    for j = 1:sample
        [tout, Nout] = uniformExact(tstop, N0, forces(i), kappa, k0b,eta,kc,ks);
        tmean(i) = tmean(i) + tout(end-1);
    end
    tmean(i) = tmean(i)/sample;
end
figure
plot(forces, tmean)