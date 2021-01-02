% William Ojemann
% BioE 485 Final Project
function [weights, ICMatrix] = bondNumberFinder(Aw)
% This function outputs a matrix of the initial conditions for the
% kolmogorov forward equations. It calculates the distribution of initial
% conditions based on the 
%   Aw - float - bond adhesion probability in the form of a decimel.

% The following is taken from Snook et al.
% For high density Aw is 0.58
% For low density Aw is 0.25
%% Generating Figures
% format long
% clear all; close all; clc;
%  k = linspace(1,10,10);
%  Pk = @(k,Aw) (1-Aw)*(-log(1-Aw)).^k./factorial(k);
%  figure
%  hold on
%  plot(k,Pk(k,0.95),'o--',k, Pk(k,0.25),'^--',k,Pk(k,0.58),'<--','LineWidth',2)
%  set(gca,'FontSize',12)
%  title('Distribution of Initial Condition Probabilities')
%  xlabel('Number of initial bonds')
%  ylabel('Probability')
%  legend('High Site Density', 'Low Site Density')
% x = linspace(1,10,10);
% Y = @(Aw) poisspdf(1:10,-log(1-Aw))*1000;
% figure
% bar(x,Y(0.58),1)
% set(gca,'FontSize',12)
% title('Distribution of Initial Bond Numbers (High Density)')
% ylabel('Number of simulations')
% xlabel('Number of initial bonds')
% figure
% bar(x,Y(0.25),1)
% set(gca,'FontSize',12)
% title('Distribution of Initial Bond Numbers (Low Density)')
% ylabel('Number of simulations')
% xlabel('Number of initial bonds')
%% Output
%Aw = 0.58;
Y = @(Aw) poisspdf(1:7,-log(1-Aw))*1000;
initials = [0 round(Y(Aw))];
ICMatrix = [zeros(1,7); diag(diag(ones(7)))];
weights = (initials(initials > 0));
%weights = weight./(sum(weight))
N = length(weights);
ICMatrix = ICMatrix(:,1:N);


% initialBonds = zeros(8,sum(initials));
% count = 2;
% spot = 1;
% for i = initials(1:7)
%     for j = 1:i
%         initialBonds(count,spot) = 1;
%         spot = spot + 1;
%     end
%     count = count + 1;
% end
save('initial_conditions.mat','initials','-mat')
end