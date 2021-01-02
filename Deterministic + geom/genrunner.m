% 5/28/20
% general script to run varied lifetime

%std deviation numbers
stdata = [5.981995393982141E-12   1.117118131404848
  11.879121014737230E-12   1.258281889845138
  18.041235762195786E-12   1.359373827229811
  29.853021931033279E-12   1.360952954674168
  39.127711940341570E-12   1.358182170871337];

% estimate=[kappa, eta, ks, kc, k01, D];
estimate = [0.00124402839199090,2.55295462419621e-10,0.111384250721051,1.56374425143498,1.18318602769146e-12,8.03092771313965e-11];

%aw is adhesion frequency .58  is HD and .25 LD
aw = .58;

%fspan is the vector of forces
fspan = stdata(:,1);
ltime = variedIClifetime(estimate,aw, fspan);
plot(stdata(:,1),ltime,stdata(:,1), stdata(:,2));
