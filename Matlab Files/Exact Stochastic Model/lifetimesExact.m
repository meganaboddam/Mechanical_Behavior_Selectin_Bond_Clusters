function [forces,tmean] = lifetimesExact(param)
%Function to calculate the mean lifetimes using the exact time step
%solution from Gillespie's 1970 work.
%   param - array - every parameter used to be fit to the objective
%   function


tstop = 200;
sample = 50;
forces = linspace(1,50)*1E-12;
tmean = zeros(1,length(forces));
for i = 1:length(forces)
    for j = 1:sample
        [tout, ~] = uniformExact(tstop, param(6), forces(i), param(1), param(2),param(3),param(4),param(5));
        tmean(i) = tmean(i) + tout(end-1);
    end
    tmean(i) = tmean(i)/sample;
end
end

