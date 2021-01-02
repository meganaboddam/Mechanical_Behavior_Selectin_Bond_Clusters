% Will Ojemann 
%BioE 485 Final Projects
% Bond Cluster Group
clearvars;close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%--High Density (cluster) data--%%%%%%%%%%%%%%%%%%
% Mean Lifetime
HDlife = [5.9519E-12  1.3129
        12.0350E-12    1.3329
        18.0890E-12    1.8626
        29.8906E-12    1.3289
        39.2560E-12    1.3382];

% Standard Deviation
HDstd= [5.981995393982141E-12   1.117118131404848
       11.879121014737230E-12   1.258281889845138
       18.041235762195786E-12   1.359373827229811
       29.853021931033279E-12   1.360952954674168
       39.127711940341570E-12   1.358182170871337];
       
%%%%%%%%%%%%%%%%%%%%%%%%%--Low Density (single) data--%%%%%%%%%%%%%%%%%%%%
% Mean Lifetime
LDlife= [6.04805000000000e-12,0.433328000000000;
         1.20022000000000e-11,0.806065000000000;
         1.80241000000000e-11,0.912796000000000;
         2.98636000000000e-11,0.799691000000000;
         3.91329000000000e-11,0.296444000000000];

%Standard Deviation
LDstd= [6.09917e-12	0.635847; 
         12.2032e-12 1.06711;
         18.1886e-12	1.16884;
         30.1231e-12	1.01431;
         39.2188e-12	0.435177];

%std deviation numbers
stdata = [5.981995393982141E-12   1.117118131404848
  11.879121014737230E-12   1.258281889845138
  18.041235762195786E-12   1.359373827229811
  29.853021931033279E-12   1.360952954674168
  39.127711940341570E-12   1.358182170871337];
%%%%%%%%%%%%%%%%%%%%%%%%%--Plotting--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check1 = 0;
while check1 == 0
    check = input('Please indicate the figure number you want to view: ');
    close all
switch check
    case 13
        % HD Mean Fit no Out
        figure
        HDlife1 = [5.9519E-12  1.3129
                12.0350E-12    1.3329
                %18.0890E-12    1.8626
                29.8906E-12    1.3289
                39.2560E-12    1.3382];
        %estimates = [0.00082089  1.5426e-10    0.097669      2.1255      3.6073  1.4126e-08];
        estimates = [-2.1704e-05   1.768e-12  -0.0085935     0.93917      6.6954  1.4914e-08];
        [time,~] = variedIClifetime2(estimates, 0.58, HDlife(:,1));
        plot(HDlife(:,1)/10^-12, time, '*', 'LineWidth', 1); hold on;
        plot(HDlife(:,1)/10^-12, HDlife(:,2),'o', 'MarkerSize', 10, 'lineWidth',1); hold off;
        ylim([0 2]);
        set(gca,'FontSize',14)
        title("Model Fit to High Density Mean Data");
        xlabel("Force (N)");
        ylabel("Bond Lifetime (s)");
        legend("Model Mean Lifetime", "Experimental Mean Lifetime");  
        %% LD Mean Fit no Out
    case 14     
        figure
        %estimates = [0.00082089  1.5426e-10    0.097669      2.1255      3.6073  1.4126e-08];
        estimates = [-2.1704e-05   1.768e-12  -0.0085935     0.93917      6.6954  1.4914e-08];
        [time1,~] = variedIClifetime2(estimates, 0.25, LDlife(:,1));
        plot(LDlife(:,1), time1, '*'); hold on;
        plot(LDlife(:,1), LDlife(:,2),'o', 'MarkerSize', 10); hold off;
        ylim([0 2]);
        set(gca,'FontSize',14)
        title('Model Fit to High Density Mean Data, Low Densit ICs');
        xlabel("Force (N)");
        ylabel("Bond Lifetime (s)");
        legend("Model Mean Lifetime", "Experimental Mean Lifetime");
        %% Bond Probabilities
    case 7
        figure
        estimates = [2.833E-3 3.73E-10 6.45E-2 2.590 6.491E-14 1.319E-09];
        [T,Y] = ode45(@BondProbability10,[0 8], [0 0 0 0 0 1 0 0], [],[estimates 18.1886e-12]);
        plot(T,Y(:,1),'--','LineWidth',2)
        hold on
        plot(T,Y(:,2:8),'LineWidth',2)
        hold off
        set(gca,'FontSize',14)
        title('Probability of Bound Receptor-Ligand Pair Cluster')
        xlabel('Time (s)')
        ylabel('Probability')
        ylim([-0.2 1])
        legend('0 bonds', '1 bond','2 bonds','3 bonds','4 bonds','5 bonds','6 bonds','7 bonds')
        %% IC Probabilities
    case 2
        k = linspace(1,10,10);
         Pk = @(k,Aw) (1-Aw)*(-log(1-Aw)).^k./factorial(k);
        figure
        semilogy(k, Pk(k,0.58),'o--','LineWidth',2)
         hold on
         semilogy(k, Pk(k,0.25),'o--','LineWidth',2)
         hold off
         set(gca,'FontSize',14)
         title('Distribution of Initial Condition Probabilities')
         xlabel('Number of initial bonds')
         ylabel('log(Probability)')
         legend('High Site Density', 'Low Site Density')
        %% ICs High Density
    case 3
        figure
        x = linspace(1,10,10);
        Y = @(Aw) poisspdf(1:10,-log(1-Aw))*1000;
        bar(x,round(Y(0.58)),1)
        set(gca,'FontSize',14)
        title('Distribution of Initial Bond Numbers for High Density ICs')
        ylabel('Number of simulations')
        xlabel('Number of initial bonds')
        %% ICs Low Density
    case 4
        figure
        x = linspace(1,10,10);
        Y = @(Aw) poisspdf(1:10,-log(1-Aw))*1000;
        bar(x,Y(0.25),1)
        set(gca,'FontSize',14)
        title('Distribution of Initial Bond Numbers Low Density ICs')
        ylabel('Number of simulations')
        xlabel('Number of initial bonds')
        %% Single Bond Lifetime
    case 5
        figure
        estimates = [2.833E-3 3.73E-10 6.45E-2 2.590 6.491E-14 1.319E-09];
        [means,~] = variedIClifetime2(estimates,0.25,(0:1:40)*10^-12);
        plot((0:1:40),means,'LineWidth',2)
        set(gca,'FontSize',14)
        ylim([0 2])
        title('Model of Single Catch Bond')
        ylabel('Lifetime (s)')
        xlabel('Force (pN)')
        %% Single Bond Insensitive
    case 6
        figure
        D = 0; 
        eta = 0;
        kc = 2.59;
        ks = 6.45e-2;
        k01 = 6.491e-14;
        kappa = 2.833e-3;
        guesses2 = [kappa, eta, ks, kc, k01, D];
        f = linspace(0, 40e-12, 50);
        [time2,~] = variedIClifetime2(guesses2, 0.25, f);
        plot(f/10^-12, time2, 'LineWidth', 2);
        ylim([0 2]);
        set(gca,'FontSize',14)
        title("Model of Force Insensitive Single Catch Bond");
        xlabel("Force (pN)");
        ylabel("Lifetime (s)");
        %% High Density Lifetime
    case 11
        figure
        estimates = [0.000249063937130034,2.46736175110646e-10,0.129473114091613,1.69761749456730,0.00185002337503608,.009e-6];
        fspan2=linspace(0,43e-12,40);
        fs = fspan2./1e-12;
        fss = HDstd(:,1)./1e-12;
        fsss = HDlife(:,1)./1e-12;
        [ltime1,~] = variedIClifetime2(estimates,.58, HDlife(:,1));
        plot(fss,HDstd(:,2),'^',fsss,HDlife(:,2),'o',fsss,ltime1,'k*')
        ylim([0 2])
        xlabel 'Force (pN)'
        ylabel 'Lifetime (s)'
        title 'High Density Cluster Bond Lifetime'
        legend ('Experimental Standard Deviation', 'Experimental Mean Lifetime', 'Model Mean Lifetime')
        set(gca,'FontSize',14)
        %% Low Density Lifetime
    case 12
        figure
        estimates = [0.000249063937130034,2.46736175110646e-10,0.129473114091613,1.69761749456730,0.00185002337503608,.009e-6];
        fspan2=linspace(0,43e-12,50);
        fs = fspan2./1e-12;
        fss = LDstd(:,1)./1e-12;
        fsss = LDlife(:,1)./1e-12;
        [ltime2,~] = variedIClifetime2(estimates,.25,LDlife(:,1));
        plot(fss,LDstd(:,2),'^',fss,LDlife(:,2),'o',fsss,ltime2,'k*')
        ylim([0 2])
        xlabel 'Force (pN)'
        ylabel 'Lifetime (s)'
        title 'Low Density Bond Lifetime'
        legend ('Experimental Standard Deviation', 'Experimental Mean lifetime', 'Model Mean Lifetime')
        set(gca,'FontSize',14)
        %% Residuals
    case 15
        figure
        fspan2=linspace(0,43e-12,40);
        fs = fspan2./1e-12;
        fss = HDstd(:,1)./1e-12;
        fsss = HDlife(:,1)./1e-12;
        estimates = [0.000249063937130034,2.46736175110646e-10,0.129473114091613,1.69761749456730,0.00185002337503608,.009e-6];
        [ltime1,~] = variedIClifetime2(estimates,.58, HDlife(:,1));
        [ltime2,~] = variedIClifetime2(estimates,.25,LDlife(:,1));
        r = (HDlife(:,2)' - ltime1)./ltime1;
        r2 = (LDlife(:,2)' - ltime2)./ltime2;
        plot(fsss, r, 'o-')
        hold on
        plot(fss, r2, '*-')
        title 'Weighted Residuals' 
        xlabel 'Force (pN)'
        ylabel 'Weighted Residual'
        legend( 'High Density' ,'Low Density')
        set(gca,'FontSize',14)
        ylim([-1 1])
        %% Spring Constant Sensitivity
    case 16
        figure
        param = [.001 .01 .1 1 10 100];
        keys = {'ro-','gd--','bs-.','c^','k>','m<'};
        estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.009e-6];
        hold on
        for i = 1:length(param)
           estimates(1) = param(i);
           [testtime,~] = variedIClifetime2(estimates,0.58,(1:1:40)*10^-12);
           plot((1:1:40),testtime,keys{i},'LineWidth', 2)
        end
        hold off
        ylim([0 2])
        set(gca,'FontSize',14)
        legend('\kappa = 0.001 N/m','\kappa = 0.01','\kappa = 0.1','\kappa = 1','\kappa = 10','\kappa = 100')
        title('Sensitivity of Model to \kappa')
        xlabel('Force (pN)')
        ylabel('Mean Lifetime (s)')
        %% K0r Sensitivity
    case 17
        figure
        param = [.001 .01 .1 1 100 1000];
        keys = {'ro-','gd--','bs-.','c^-','k>--','m<-.'};
        estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.009e-6];
        hold on
        for i = 1:length(param)
           estimates(5) = param(i);
           [testtime,~] = variedIClifetime2(estimates,0.58,(1:1:40)*10^-12);
           plot(((1:1:40)),testtime,keys{i},'LineWidth', 2)
        end
        hold off
        %ylim([0 2])
        set(gca,'FontSize',14)
        legend('k0r = 0.001 1/s','k0r = 0.01','k0r = 0.1','k0r = 10','k0r = 100','k0r = 1000')
        title('Sensitivity of Model to K0r')
        xlabel('Force (pN)')
        ylabel('Mean Lifetime (s)')
        %% Xi Sensitivity
    case 10
        figure
        param = [1 10 20 30 40 50]*10^-11;
        keys = {'ro-','gd--','bs-.','c^-','k>--','m<-.'};
        estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.009e-6];
        hold on
        for i = 1:length(param)
           estimates(2) = param(i);
           [testtime,~] = variedIClifetime2(estimates,0.58,(1:1:40)*10^-12);
           plot((1:1:40),testtime,keys{i},'LineWidth', 2)
        end
        hold off
        ylim([0 2])
        set(gca,'FontSize',14)
        legend('\xi = 1E-11 m','\xi = 10E-11','\xi = 20E-11','\xi = 30E-11','\xi = 40E-11','\xi = 50E-11')
        title('Sensitivity of Model to \xi')
        xlabel('Force (pN)')
        ylabel('Mean Lifetime (s)')
        %% Ks Sensitivity
    case 8
        figure
        param = [.001 .01 .1 1 10 100];
        keys = {'ro-','gd--','bs-.','c^-','k>--','m<-.'};
        estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.009e-6];
        hold on
        for i = 1:length(param)
           estimates(3) = param(i);
           [testtime,~] = variedIClifetime2(estimates,0.58,(1:1:40)*10^-12);
           plot((1:1:40),testtime,keys{i},'LineWidth', 2)
        end
        hold off
        set(gca,'FontSize',14)
        legend('k_s = 0.001 (1/s)','k_s = 0.01','k_s = 0.1','k_s = 1','k_s = 10','k_s = 100')
        title('Sensitivity of Model to k_s')
        xlabel('Force (pN)')
        ylabel('Mean Lifetime (s)')
        %% Kc Sensitivity
    case 9
        figure
        param = [.001 .01 .1 1 10 100];
        keys = {'ro-','gd--','bs-.','c^-','k>--','m<-.'};
        estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.009e-6];
        hold on
        for i = 1:length(param)
           estimates(4) = param(i);
           [testtime,~] = variedIClifetime2(estimates,0.58,(1:1:40)*10^-12);
           plot((1:1:40),testtime,keys{i},'LineWidth', 2)
        end
        hold off
        set(gca,'FontSize',14)
        legend('k_c = 0.001 (1/s)','k_c = 0.01','k_c = 0.1','k_c = 1','k_c = 10','k_c = 100')
        title('Sensitivity of Model to k_c')
        xlabel('Force (pN)')
        ylabel('Mean Lifetime (s)')
        % %% D Sensitivity
        % figure
        % param = [.001 .01 .1 1 10 100]*10^-8;
        % keys = {'ro-','gd--','bs-.','c^-','k>--','m<-.'};
        % estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.006e-6];
        % hold on
        % for i = 1:length(param)
        %    estimates(6) = param(i);
        %    [testtime,~] = variedIClifetime2(estimates,0.58,(1:1:40)*10^-12);
        %    plot((1:1:40)*10^-12,testtime,keys{i},'LineWidth', 2)
        % end
        % hold off
        % ylim([0 2])
        % set(gca,'FontSize',14)
        % legend('D = .001E-8 m','D = .01E-8','D = .1E-8','D = 1E-8','D = 10E-8','D = 100E-8')
        % title('Sensitivity of Model to D')
        % xlabel('Force (N)')
        % ylabel('Mean Lifetime (s)')
        %% kon force sensitivity
    case 18
        estimates = [0.000249063937130034, 2.46736175110646e-10, 0.129473114091613, 1.69761749456730, 0.00185002337503608,0.009e-6];
        figure
        rates = sun(estimates(1),estimates(5),estimates(6),[1 2 3 4 5 6 7], HDlife(:,1));
        plot(data(:,1)*10^12, rates,'o-','LineWidth',2)
        set(gca,'FontSize',14)
        legend('1 bond', '2 bonds', '3 bonds', '4 bonds', '5 bonds', '6 bonds', '7 bonds')
        title('Sensitivity of Rebinding Rate to Force')
        xlabel('Force (pN)')
        ylabel('Rebinding Rate (1/s)')
    otherwise
        disp('Sorry, not a figure number :( try again!')
end % Switch
mike = input('Do you want to plot another figure (y/n)? ', 's');
if mike == 'n'
    check1 = 1;
elseif mike == 'y'
    check1 = 0;
else
    disp("Not a valid answer, you're going again!")
end % if
end % While