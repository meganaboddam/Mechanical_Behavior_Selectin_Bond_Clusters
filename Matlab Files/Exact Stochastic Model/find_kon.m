function [kon] = find_kon(kappa,n,force)
%Finds the onrate of a catch bond in a cluster
%   kappa is spring constant(units:N/nm),
%   n is number of ligands bonded to receptors,
%   force is applied constantly, kon calculation taken from Whitfield


k01 = 4; %zero distance on rate, greater than 3 based on whitfield 
kT = 298*1.38*10^-23; %bolz constant*room temp(Joules)
d = force/(n*kappa);
xrms = sqrt(kT/kappa);
kon = k01*erfc(d/(sqrt(2)*xrms));
end