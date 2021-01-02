%Bond Cluster Group 485 
%Objective function for cluster bonds
function J = objj(params)
%   params - the six parameters that define the different on and off
%   equations in the format [kappa, eta, ks, kc, k01, D]

%------------------------------------------------------------------------
%                           DATA FOR COMPARISON AGAINST
%     
%
%                           Data format - [ force, lifetime]
%
%                         --High Density (Cluster) data -- 
% Mean Lifetime
HDlife = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        18.0890E-12    1.8626
        29.8906E-12    1.3289
        39.2560E-12    1.3382];

% Standard Deviation
HDstd= [5.981995393982141E-12   1.117118131404848
       11.879121014737230E-12   1.258281889845138
       18.041235762195786E-12   1.359373827229811
       29.853021931033279E-12   1.360952954674168
       39.127711940341570E-12   1.358182170871337];
       
%                          --Low Density (single) data --
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


 

% to kepp track when parameter fitting
disp(params)

% gives the lifetime for the low density (single) model
%Change between LDlife and LDstd to fit differently
[time2,~] = variedIClifetime2(params,.25, LDlife(:,1));
%[nan,std2] = variedIClifetime2(params,.25, LDstd(:,1));

% gives the lifetime for the high density 
% change beteween HDlife and HDstd to fit differently
[time,~] = variedIClifetime2(params,.58,HDlife(:,1));
%[nan,std] = variedIClifetime2(params,.58,HDstd(:,1));

%plots the current fit of the data
plot(HDlife(:,1),time,'*',HDlife(:,1),HDlife(:,2),'o',LDlife(:,1),LDlife(:,2),'d'); drawnow


%objective function
%change obj function to whatever take into account both systems or just one

J = (sum((HDlife(:,2)'-time).^2))+(sum((LDlife(:,2)'-time2).^2));%+(sum((HDstd(:,2)'-std).^2))+(sum((LDstd(:,2)'-std2).^2));
end
    

