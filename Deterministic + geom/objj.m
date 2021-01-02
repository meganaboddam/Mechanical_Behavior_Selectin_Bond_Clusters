function J = objj(params)
% original lifetime data
data = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        18.0890E-12    1.8626
        29.8906E-12    1.3289
        39.2560E-12    1.3382];
% standard deviation data
%   data_std = [5.981995393982141E-12   1.117118131404848
%   11.879121014737230E-12   1.258281889845138
%   18.041235762195786E-12   1.359373827229811
%   29.853021931033279E-12   1.360952954674168
%   39.127711940341570E-12   1.358182170871337];
%    data2 = [
%   18.041235762195786E-12   1.359373827229811
%   29.853021931033279E-12   1.360952954674168
%   39.127711940341570E-12   1.358182170871337];
%low density single bond data
%dataLD = importdata('lowDensity_data.dat');
%dataLD(:,1) = dataLD(:,1).*(1E-12);
%disp(params)

% gives the lifetime for the singlebond system
%time2 = variedIClifetime(params,.25, dataLD(:,1))';

% kc = 10;
% k01 = 1e-12; % rate constant from sun et al
% D = .01e-8;
% kappa = 0.00131211798645895;
% eta = 1.51813939524036e-10;
% ks = 0.255082466895303;


%params = [kappa, eta, ks, kc, k01, D];
% gives the lifetime for the varied initial condition 7 state model
% 0.58 is adhesion frequency



[time, std] = variedIClifetime(abs(params),.58,data(:,1));


%plots the current fit of the data
plot(data(:,1),time,'*',data(:,1),data(:,2),'o'); drawnow
%objective function
%change obj function to whatever you are running for
%

wres = sum(((data(:,2)-time).^2) ./ std);
J = sum(wres); 


%+(sum((dataLD(:,2)-time2).^2));
end
    

