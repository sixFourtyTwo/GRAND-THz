
figure(1)

% Plot for ZF channel equalization
% plot1 = semilogy(p.SNRdB_list, sum(ber_turbo_ZF_hard_Unc,2)/p.trials, 'b--x');
% hold on

% plot2 = semilogy(p.SNRdB_list, sum(ber_turbo_ZF_hard,2)/p.trials, 'm--*');
% hold on
% % 
% % plot3 = semilogy(p.SNRdB_list, sum(ber_turbo_ZF_soft_Unc,2)/p.trials, 'k--s');
% % hold on
% 
% plot4 = semilogy(p.SNRdB_list, sum(ber_turbo_ZF_soft,2)/p.trials, 'g--o');
% hold on

% Plot for ML channel equalization
% plot5 = semilogy(p.SNRdB_list, sum(ber_turbo_ML_hard_Unc,2)/p.trials, 'r--^');
% hold on

plot6 = semilogy(p.SNRdB_list, sum(ber_turbo_ML_hard,2)/p.trials, 'c--d');
hold on

% plot7 = semilogy(p.SNRdB_list, sum(ber_turbo_ML_soft_Unc,2)/p.trials, 'y--^');
% hold on

plot8 = semilogy(p.SNRdB_list, sum(ber_turbo_ML_soft,2)/p.trials, 'b--d');
hold on

% Plot for SP channel equalization
% plot9 = semilogy(p.SNRdB_list, sum(ber_turbo_SP_hard_Unc,2)/p.trials, 'r--^');
% hold on

plot10 = semilogy(p.SNRdB_list, sum(ber_turbo_SP_hard,2)/p.trials, 'r--^');
hold on

% plot11 = semilogy(p.SNRdB_list, sum(ber_turbo_SP_soft_Unc,2)/p.trials, 'y--^');
% hold on

plot12 = semilogy(p.SNRdB_list, sum(ber_turbo_SP_soft,2)/p.trials, 'y--^');
hold on

% Set labels and legends
xlabel('SNR (dB)');
ylabel('BER');
% legend([plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12], ...
%     'ZF Hard (Unc)', 'ZF Hard', 'ZF Soft (Unc)', 'ZF Soft', 'ML Hard (Unc)', 'ML Hard', 'ML Soft (Unc)', 'ML Soft', 'SP Hard (Unc)', 'SP Hard', 'SP Soft (Unc)', 'SP Soft');
% grid on

legend([ plot6,  plot8,  plot10, plot12], ...
        'ML Hard', 'ML Soft', 'SP Hard',  'SP Soft');


