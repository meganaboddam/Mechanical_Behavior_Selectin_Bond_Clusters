function dPdt = BondProbability7Geom(t,p, param)
% force,kOff_without_force,kappa,kT,dxOff,n
%BondProbability(p)
%BondProbability(p,force,kOff_initial,kappa,n)
%function to solve the ODEs from the chemical Model
  % n = [1 2 3 4 5 6 7 8 9 10];
  n = [1 2 3 4 5 6 7];

 
% parameters for kon and koff along with force at the end
 kappa = param(1);
 eta = param(2);
 ks = param(3);
 kc = param(4);
%  k10 = 1e-12; % rate constant from sun et al
% D = .01e-6;
 k10 = param(5);
 D = param(6);
 force = param(7);
 
 
 %kappa = 10;
% these parameters are estimates from Novikova et al.

% produces forces for middle and outside geometries
x = 1.0E-10;
[fMiddle,fOutside] = forceDistributer(x,force,n,kappa);
% produces vector of kOff corresponding to bound states
kOff = findkoffGeom(fMiddle,eta,ks,kc, n);
kOffOut = findkoffGeom(fOutside,eta,ks,kc, n);
% produces vector of kon corresponding to bound states
kon = sunGeom2(fMiddle, k10, D, kappa, n);
konOut = sunGeom2(fOutside, k10, D, kappa, n);

    dPdt = zeros(8,1);%creates an empty matrix where the solution will go
    % can be altered to run for 10 bound states
    dPdt(1)= p(2) * kOff(1); % 0 bonds
    
    dPdt(2)=  p(3) .* (kOff(2) + kOffOut(2)) - p(2) .* kon(1) - p(2) .* kOff(1); %1 bonds 
    
    dPdt(3)= -1 * p(3) .* (kOff(2) + kOffOut(2)) - 1 * p(3) .* (kon(2) + konOut(2)) + p(2) .* kon(1) + p(4) .* (kOff(3) + 2 * kOffOut(3)); % 2 bonds
    
    dPdt(4)= -1 * p(4) .* (kOff(3) + 2 * kOffOut(3)) - 1 * p(4) .* (kon(3) + 2 * konOut(3)) + p(3) .* (kon(2) + konOut(2)) + p(5) .* (kOff(4) + 3 * kOffOut(4)); %3 bonds
    
    dPdt(5)= -1 * p(5) .* (kOff(4) + 3 * kOffOut(4)) - 1 * p(5) .* (kon(4) + 3 * konOut(4)) + p(4) .* (kon(3) + 2 * konOut(3)) + p(6) .* (kOff(5) + 4 * kOffOut(5)); %4 bonds
    
    dPdt(6)= -1 * p(6) .* (kOff(5) + 4 * kOffOut(5)) - 1 * p(6) .* (kon(5) + 4 * konOut(5)) + p(5) .* (kon(4) + 3 * kon(4)) + p(7) .* (kOff(6) + 5 * kOffOut(6)); %5 bonds
    
    dPdt(7)= -1 * p(7) .* (kOff(6) + 5 * kOffOut(6)) - 1 * p(7) .* (kon(6) + 5 * konOut(6)) + p(6) .* (kon(5) + 4 * konOut(5)) + p(8) .* (kOff(7) + 6 * kOffOut(6)); %6 bonds
    
    dPdt(8)= -1 * p(8) .* (kOff(7) + 6 * kOffOut(7)) + p(7) .* (kon(6) + 5 * konOut(6)); %7 bonds
    
    
   
    
    
    
    %-7*p(8).*kon(7)+8.*p(9).*kOff(8); %7 bonds
    
    
    
    
    
    
%     dPdt(9)= -8 * p(9).* kOff(8) -8*p(9).*kon(8)+ 7*p(8) .* kon(7)+ 9.*p(10).*kOff(9); %8 bonds
%     dPdt(10)= -9 * p(10).* kOff(9) -9*p(10).*kon(9)+ 8*p(9) .* kon(8)+ 10.*p(11).*kOff(10); %9 bonds
%     dPdt(11)= -10 * p(11).* kOff(10)+ 9*p(10) .* kon(10); %10 bonds
end
