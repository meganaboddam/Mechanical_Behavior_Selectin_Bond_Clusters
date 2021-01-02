
%Bioen 458
function [kon] = sun(kappa,k01,n,D,force)
%Finds the onrate of a catch bond in a cluster
%   kappa -  is spring constant(units:N/nm)
%   k01 - rate constant for adhesion rate   
%   n - vector  that consists of bond numbers [1 2 3 4 5 6 7]
%   force - the force that is pulling the bead
%   D - characteristic distance of the ligand receptor pair, should be on
%   the near .001e-6

%k01 = 1e-12; % rate constant from sun et al
%D = .008e-6;


d = (force./n)./kappa;


kon = k01.*exp(-((d./D).^2));
end