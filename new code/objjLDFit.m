% Adapted by Megana Boddam from Frank's Code
% Bond Cluster Group 485 
% Objective function for cluster bonds
function J = objjLDFit(params)
%   params - the six parameters that define the different on and off
%   equations in the format [kappa, eta, ks, kc, k01, D]

%------------------------------------------------------------------------
%                           DATA FOR COMPARISON AGAINST
%
%                           Data format - [ force, lifetime]
%   
%                           --Low Density (single) data --
% Mean Lifetime
LDlife= [6.04805000000000e-12,0.433328000000000;
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
%
%-----------------------------------------------------------------------     
     
% to kepp track when parameter fitting
disp(params)

% gives the lifetime for the low density clusters
% change between LDlife and LDstd to fit differently
[time,~] = variedIClifetime2(params,.25, LDlife(:,1));
%[~,std2] = variedIClifetime2(params,.25, LDstd(:,1));

% objective function
J = sum(sum(((LDlife(:,2)'-time).^2)./(LDstd(:,2))));

end
    

