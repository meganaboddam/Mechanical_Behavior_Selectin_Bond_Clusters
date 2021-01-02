% William Ojemann (David, Frank, Jasmeet, Megana)
% BioE 485 Final Project
% 5/21/2020
function [koff] = findkoff(N,f,eta,ks,kc)
% This function finds the koff rate of unbinding for a catch bond cluster
% model based on the number of closed bonds N and force f.
% This work comes from Novikova and Storm 2013
%   N - int - number of closed bonds
%   f - int - force applied to the cluster
%   k0b - float - base unbinding rate (to be fit)
%   eta - float - parameter for unbinding lengths (to be fit)
%   pc - float - unified critical force for catch bond pathway (to be fit)
%   ps - float - unified critical force for slip bond pathway (to be fit)
%   fstar - float - unifying parameter equal to kB*T/eta where eta is
%   characteristic bond length (to be fit)

kB = 1.38*10^-23; % Boltzman Constant
T = 298; % Temperature (K)
f = f./N;
koff = (kc*exp(-(f.*eta/(kB*T)))) + ks*exp(f.*eta/(kB*T));
end

