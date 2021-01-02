%CatchBondClusterFive
clf
clc
clear

data1 = [5.9519E-12, 12.0350E-12, 18.0890E-12, 29.8906E-12, 39.2560E-12];
data2 = [1.3129, 1.3329, 1.8626, 1.3289, 1.3382];


plot (data1,data2, '*')
hold on 
forces = [6E-12, 12E-12, 18E-12, 30E-12, 39E-12];
for i = 1:length(forces)
n = 5;
numberBonds = [1, 2, 3, 4, 5];
%totalForce = 10000E-12; %Pn
forcePerBond = forces(i)./numberBonds;
kOff_without_force = 1.3E-7;
kappa = 2.4E-6;
kT = 4.11 * 10 ^ (-21);
dxOff = .3 * 10 ^ (-9);

k0b = 1;
eta = 1E-8;
kc = 1E-7;
ks = 1E-5;

%kOffArray = find_koffWill(kT,forcePerBond,kOff_without_force,dxOff);
kOffArray = find_koffWill(numberBonds,forces(i),k0b,eta,kc,ks);
%kOffArray = find_kOff_single_bond(kT,forcePerBond,kOff_without_force,dxOff);
kOnArray = sun(kappa,forcePerBond);

P5i = 1; % intially there are 5 bonds formed
P4i = 0; % all bond are formed
P3i = 0; % all bond are formed
P2i = 0; % all bond are formed
P1i = 0; % all bond are formed
P0i = 0; % all bonds are formed initially
%IC = [P0i, P1i, P2i, P3i, P4i, P5i];
%TSPAN = 0:10E3:10E9;
IC = [0,1,0,0,0,0]
TSPAN = 0:0.1:10;

[Tout,Pout] = ode23s(@(t,p)BondProbabilityFive(p,kOffArray, kOnArray),TSPAN,IC);
Probability = Tout.*Pout;

probBound = 1 - Pout(:,1);
differprobBound = gradient(probBound)./gradient(Tout);
timet = Tout.* differprobBound;
lifetime(i) = -trapz(0.1, timet);

end

plot(forces, lifetime, 'o')
ylim ([0 2])
plot(Tout,Pout(:,1),'r')
hold on 
plot(Tout,Pout(:,2),'b')
plot(Tout,Pout(:,3),'g')
plot(Tout,Pout(:,4),'m')
plot(Tout,Pout(:,5),'k')
plot(Tout,Pout(:,6),'c')
legend("No bound", "One bound", "Two bound", "Three bound", "Four Bound", "Five bound");
xlabel('Time')
ylabel("Probability")
title("Probability of Bound Catch Bonds in a Cluster")
 
 
% % deriv = diff(Pout(:,5))./diff(Tout);
% % integ = trapz(deriv, TSPAN(1,1000000));
% 
% lamda = 5/10;
% pd = makedist('Poisson','lambda',lamda);
% pdf_normal = pdf(pd,Pout);
% 
% %forces = 
% variable = Tout .* (1-pdf_normal(:,1));
% lifetime = trapz(0.1,variable);
% %plot(Tout,lifetime)