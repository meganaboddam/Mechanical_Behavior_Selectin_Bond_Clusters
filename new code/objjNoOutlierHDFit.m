% Adapted by Megana Boddam from Frank's Code
% Bond Cluster Group 485 
% Objective function cluster bonds
function J = objjNoOutlierHDFit(params)
%   params - the six parameters that define the different on and off
%   equations in the format [kappa, eta, ks, kc, k01, D]

%------------------------------------------------------------------------
%                           DATA FOR COMPARISON AGAINST
%
%                           Data format - [ force, lifetime]
%
%                         --High Density (Cluster) data -- 
%       
% Data after removing the outlier
HDlife = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        29.8906E-12    1.3289
        39.2560E-12    1.3382];
HDstd= [5.981995393982141E-12   1.117118131404848
           11.879121014737230E-12   1.258281889845138
           29.853021931033279E-12   1.360952954674168
           39.127711940341570E-12   1.358182170871337]; 
%       
%-------------------------------------------------------------------------
       
% to keep track while parameter fitting
disp(params)

% gives the lifetime for the high density 
% change beteween HDlife and HDstd to fit differently
[time,~] = variedIClifetime2(params,.58,HDlife(:,1));
%[~,std] = variedIClifetime2(params,.58,HDstd(:,1));

% objective function
J = sum(sum(((HDlife(:,2)'-time).^2)./(HDstd(:,2))));

end
    

