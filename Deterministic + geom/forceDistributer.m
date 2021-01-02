function [fMiddle,fOutside] = forceDistributer(x,F,N,k)
%forceDistributer takes in the parameters above and returns the force on
%the individual bond in the middle of the geometry and each individual bond
%on the edge of the geometry. In the case of N = 1 the force on the single
%bond is fMiddle.
%   x - float - maximum radius in which bonds can form to the magnetic bead
%   F - float - force distributed among all of the bonds
%   N - int vector - number of bonds
r = 1.5E-9; % radius of the magnetic bead
fOutside = F./N + x*tan(0.5*asin(x/r))*k;
fMiddle = F - fOutside .* (N-1);

% fMiddle = F./N;
% fOutside = F./N;

end