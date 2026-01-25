%% Given a codebook calculate its achievable rate

%Codebooks are consider to be already defined, this example can be run
%after running Codes_T4_M2 to obtaine Cg2, Cg4, Cg8 and Cg16


N = 2;

niter = 10^6; %Number of iteration for the Monte Carlo Simulation
snr = -6:2:20;  %snr in dBs;
snrn = length(snr); %Number of snr points
rho = 10.^(snr./10); %snr values in natural units

Ach_Rate2 = zeros(1,snrn);
Ach_Rate4 = zeros(1,snrn);
Ach_Rate8 = zeros(1,snrn);
Ach_Rate16 = zeros(1,snrn);
Ach_Rate16ub = zeros(1,snrn);
capacity = zeros(1,snrn);

for n=1:snrn
    Ach_Rate2(n)= achievable_rate_noncoherent(Cg2,N,rho(n),niter);
    Ach_Rate4(n)= achievable_rate_noncoherent(Cg4,N,rho(n),niter);
    Ach_Rate8(n)= achievable_rate_noncoherent(Cg8,N,rho(n),niter);
    Ach_Rate16(n)= achievable_rate_noncoherent(Cg16,N,rho(n),niter);
    Ach_Rate16ub(n)= achievable_rate_noncoherent(Cg16ub,N,rho(n),niter);
end
%% Plot results

fs = 18;
lw = 1.5;
ms = 8;

figure()
plot(snr,Ach_Rate2,'m--v','MarkerSize',ms,'LineWidth',lw)
hold on
plot(snr,Ach_Rate4,'r--v','MarkerSize',ms,'LineWidth',lw)
hold on
plot(snr,Ach_Rate8,'k--v','MarkerSize',ms,'LineWidth',lw)
hold on
plot(snr,Ach_Rate16,'b--v','MarkerSize',ms,'LineWidth',lw)
hold on
plot(snr,Ach_Rate16ub,'b--*','MarkerSize',ms,'LineWidth',lw)

grid minor
ylim([0 1.5]);
xlim([min(snr) max(snr)]);
xticks([-6 0 10 20])
set(gca,'FontSize',fs)
set(gca,'TickLabelInterpreter','latex')
title('$T=4, M=2, N=2$','Interpreter', 'latex');
legend({'$L=2$','$L=4$','$L=8$','$L=16$','$L=16$ (UB)'},'Interpreter', 'latex','Location','best','NumColumns', 2);
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Achievable Rate (bpz/Hz)', 'Interpreter', 'latex');

%Save the figure 
% print(gcf, '-depsc', 'rate.eps');

