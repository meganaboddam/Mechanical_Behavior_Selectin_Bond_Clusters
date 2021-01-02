function dPdt = BondProbability(p,force,kOff_without_force,kappa,kT,dxOff,n)
%BondProbability(p)
%BondProbability(p,force,kOff_initial,kappa,n)
%function to solve the ODEs from the chemical Model
% n = 2;
% force = 150E-12;
% kOff_without_force = 1.3E-9;
% kappa = 2.4E-6;
% kT = 4.11 * 10 ^ (-21);
% dxOff = .3 * 10 ^ (-9); 

kOff = find_kOff_single_bond(kT,force,kOff_without_force,dxOff)
kon = find_kon(kappa,n,force);

    dPdt = zeros(3,1); %creates an empty matrix where the solution will go
    
    dPdt(1) = -2 * p(1) * kOff + p(2) * kon; %2 bonds
    dPdt(2)=   2 * p(1) * kOff - p(2) * kon - p(2) * kOff; %1 bonds 
    dPdt(3)= p(2) * kOff; %3 no bonds
    
end

