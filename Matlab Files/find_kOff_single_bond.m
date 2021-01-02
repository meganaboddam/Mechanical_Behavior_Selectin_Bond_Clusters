function kOff_single_bond = find_kOff_single_bond(kT,force,kOff_without_force,dxOff)

%kT = 4.11 * 10 ^ (-21);
%dxOff = .3 * 10 ^ (-9); %bond distance

% We can find k10 (unbinding rate) from Whitfield 2010 Eq. 1
kOff_single_bond_equation = @(f) kOff_without_force .* exp(dxOff*f/kT);

kOff_single_bond = kOff_single_bond_equation(force);
end

%%
% force have to specify what force it is 
% to get the force per bond == force/n
% N = sum of n
% have to calulate a different Koff and Kon depending on the force of each
% a each bond puting each kon or koff in the respective p(i) etc 
