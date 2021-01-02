function [kon] = sunGeom2(F, k01, D, kappa, N)
%Finds the onrate of a catch bond in a cluster
%   kappa is spring constant(units:N/nm),
%   n is number of ligands bonded to receptors,
%   forc3456ewe is applied constantly, kon calculation taken from sun
% n is a vector of the bond numbers [1 2 3 4 5]
%k01 = 1e-12; % rate constant from sun et al
%D = .01e-8;
% d = 1e-12:1e-10:1e-8;
d = (F ./ N) ./ (kappa);

    kon = (k01*exp(-((d./D).^2)));

end