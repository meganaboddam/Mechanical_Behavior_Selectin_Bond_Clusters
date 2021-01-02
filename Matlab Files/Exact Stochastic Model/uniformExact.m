% William Ojemann (David, Frank, Jasmeet, Megana)
% BioE 485 Final Project
% 5/21/2020
% Adapted from code produced by Harley Day:
% https://www.mathworks.com/matlabcentral/answers
% /224604-implementing-gillespie-s-algorithm
function [tout, Nout] = uniformExact(tstop, N0, f, kappa, k0b,eta,kc,ks)
% This is an exact stochastic solution to the differential equation
% describing a network of N0 independent catch bonds.
%   tstop - int - time scale for the solution
%   N0 - int - initial/total number of catch bonds
%   f - constant force for the simulation

n0 = 0; % This is the initial number of open bonds that can reclose
        % equal to N0-N

% p is normalized force with constant f* = kBT/eta
% pc is the 
% k0 is zero distance on rate for rebinding
% kappa is spring constant
% d is f/(N*kappa)
% xrms = sqrt(kB*T/kappa)
%------------------------
% equation         | rate
%------------------|-----
% N -> n           | koff = (exp(p/N-pc)+exp(p/N-ps))
% n -> N           | kon = k0*erfc(d/(sqrt(2)*xrms))
%------------------------
N = [N0; n0];
t = 0;
tout = zeros(1);
Nout = N;
iter = 0;

while t < tstop
    % getting the off rate based on the current values of N
    koff = find_koff(N(1),f,k0b,eta,kc,ks);
    % getting the on rate based on the current value of N
    kon = find_kon(kappa, N(1), f);
    % building rate array
    c = [koff, kon];
    % building species array
    h = N;
    % building reaction pathway array
    a = h'.*c;
    % generating random numbers for stochasticity
    r = rand(2,1);
    % creating random tau step for a random change in time
    tau = -log(r(1))/sum(a);
    t = t + tau;
    % creating the mechanism of choosing reaction pathway
    mu = sum(r(2)*sum(a) <= cumsum(a));
    % updating population values
    switch mu
        case 2
            N(1) = N(1) - 1;
            N(2) = N(2) + 1;
        case 1
            N(1) = N(1) + 1;
            N(2) = N(2) - 1;
    end % switch
    iter = iter + 1;
    tout = [tout t];
    Nout(:,iter+1) = N;
end % while
end % function

