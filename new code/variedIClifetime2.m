
% BioE 485 Final Project
% 5/28/2020
function [ltime,stddev] = variedIClifetime2(guess,aw, fspan)

% input
%   guess - [kappa, eta, ks, kc, k01, D] parameter guesses
%   aw - adhesion frequency: HD = .58 LD = .25
%   fspan - span of forces lifetime modeled over
 
% output
%   ltime - vector of lifetimes corresponding to forces input


%   tspan -  time span for which we are modeling over
%   ts - timestep of which we hae our time span

ts = 0.01;
tspan = 0:ts:8;

% known parameters from Snook et al. for determining min lifetime
% detectable

nu = 0.001;
r = 1.5e-6;
z = 9e-6;

%   minL -  a function that takes in force and outputs min lifetime
%   detectable by expiremental procedure

 
minL = @(fspan) (6*pi*nu*r*z)./fspan;

%   minimum - for fspan used minimum is a vector of corresponding minimum
%   lifetimes detectable

minimum = minL(fspan);

% get Initial conditions from bondNumberfinder
% inital conditions come out as a matrix where coloumns are the IC for
% different numbers of initial bonds and weights are the number of times
% each initial condition is being "run"

  [weights, ICs] =  bondNumberFinder(aw);
  
%   weights - weighting of each IC
%   ICs - vector, initial condition that gives a probability of bound states

% ---------------------------------------------------------------------------
% computing size of fspan and ICs to determine how long to run for loops to
% get all the lifetimes

  a =length(fspan);
  b = size(ICs,2);
  lifetime = zeros(b,a);
 
for j = 1:b
   %pulls first initial condition
    ic = ICs(:,j);
    %runs through all the fspan then for a single IC
for i = 1:a
   
params = [guess fspan(i)];
%calculates the prbability in the different states
[T, Y] = ode45(@BondProbability10,tspan,ic,[],params);
% translates that probability into lifetime

% calculates probability of being in all bounds state

Pbound= 1 - Y(:,1);
% calculates derivative of Pbound
dPbounddt = gradient(Pbound)./gradient(T);
% multiplies derivative by time
dpdt = T.*dPbounddt;
% Takes the integral of the derivative multiplied by t from 0 - infinitiy
life = -trapz(ts,dpdt);

% Checking if the lifetime gotten is above the minimal lifetime detecable
%If lifetime above minimum then it goes into lifetime matrix
%IF below then lifetime vector stays 0 for that spot
% if life >= minimum(i)
     lifetime(j,i)= life;
% else
% end
% plot probabilities over time 
%plot(T,Y); drawnow

end
end
% ----------------------------------------------------------------------
%   standard deviation calculation
% 

% applying the weights to the detectable lifetimes
%once weights are applied the true mean lifetime can be estimated
test = zeros(sum(weights),a);
step = zeros(1,b-1);
for i = 1:b-1
    step(i) = sum(weights(1:i))+1;
end
for i = 1:b
    for j = 1:a
        if i == 1
    test(1:weights(i),j) = lifetime(1,j).*ones(weights(i),1);
        else
            test(step(i-1):step(i-1)+weights(i)-1,j) = lifetime(i,j).*ones(weights(i),1);
        end
    end
end

% Standard deviation calculation
stand = zeros(1,a);
for i = 1:a
    lifet = test(:,i);
    lifet = lifet(lifet~=0);
    stand(1,i) = std(lifet);
end
 stddev = stand;
   
%calculates mean lifetimes over the runs for the different forces applied
ltime = mean(test);%mean(test); %sum(lifetime);
% plot( fspan, ltime)
% ylim([0 2]);
%end


