
% BioE 485 Final Project
% 5/28/2020
function [ltime, lstd]  = variedIClifetime(guess,aw, fspan)

% guess is a vector of parameters
% guess =[kappa, eta, ks, kc, k01, D];
%aw is the adhesion frequency so either .58 for HD or .25 for LD
%fspan is the span of forces put through model

% 8 second span
tspan = 0:.01:8;
 

% get Initial conditions from bondNumberfinder
% inital conditions come out as a matrix where coloumns are the IC for
% different numbers of initial bonds and weights are the number of times
% each initial condition is being "run"
  [weights, ICs] =  bondNumberFinder(aw);
  a =length(fspan);
  b = size(ICs,2);
  lifetime = zeros(b,a);
  %tins for the size of initial conditions
for j = 1:b
   %pulls first initial condition
    ic = ICs(:,j);
    %runs through all the fspan then for a single IC
for i = 1:a
   
params = [guess fspan(i)];
%calculates the prbability in the different states
[T, Y] = ode15s(@BondProbability7Geom,tspan,ic,[],params);
% translates that probability into lifetime
p5= 1 - Y(:,1);
dp5dt = gradient(p5)./gradient(T);
dpdt = T.*dp5dt;
% adding in lifetime to matrix and weighting it by initial condition.
lifetime(j,i) = trapz(.01,dpdt)*weights(j)/sum(weights);
%plot(T,Y); drawnow
end
end
%calculates mean lifetimes over the runs for the different forces applied


ltime = mean(lifetime);
lstd = std(lifetime);
% plot( fspan, ltime)
% ylim([0 2]);
%end


