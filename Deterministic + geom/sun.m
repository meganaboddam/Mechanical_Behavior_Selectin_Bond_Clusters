function [kon] = sun(kappa,k01,D,n,force)
%Finds the onrate of a catch bond in a cluster
%   kappa is spring constant(units:N/nm),
%   n is number of ligands bonded to receptors,
%   force is applied constantly, kon calculation taken from sun
% n is a vector of the bond numbers [1 2 3 4 5]
%k01 = 1e-12; % rate constant from sun et al
%D = .01e-6;
d = (force./n)./kappa;

kon = k01*exp(-((d./D).^2));
end