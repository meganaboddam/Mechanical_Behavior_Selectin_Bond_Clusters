function dPdt = BondProbabilityFive(p,kOffArray, kOnArray)
%BondProbability(p)
%BondProbability(p,force,kOff_initial,kappa,n)
%function to solve the ODEs from the chemical Model
% n = 2;
% force = 150E-12;
% kOff_without_force = 1.3E-9;
% kappa = 2.4E-6;
% kT = 4.11 * 10 ^ (-21);
% dxOff = .3 * 10 ^ (-9); 

    dPdt = zeros(6,1); %creates an empty matrix where the solution will go
    
    dPdt(6) = -5 * p(6) * kOffArray(5) +  4 * p(5) * kOnArray(4);
    dPdt(5) = 5 * p(6) * kOffArray(5) - 4 * p(5) * kOnArray(4) - 4 * p(5) * kOffArray(4) +  3 * p(4) * kOnArray(3); %1 bonds 
    dPdt(4) = 4 * p(5) * kOffArray(4) - 3 * p(4) * kOnArray(3) - 3 * p(4) * kOffArray(3) +  2 * p(3) * kOnArray(2);
    dPdt(3) = 3 * p(4) * kOffArray(3) - 2 * p(3) * kOnArray(2) - 2 * p(3) * kOffArray(2) + p(2) * kOnArray(1);
    dPdt(2) = 2 * p(3) * kOffArray(2) -     p(2) * kOnArray(1) - p(2) * kOffArray(1);
    dPdt(1) =     p(2) * kOffArray(1);
    
    
end
