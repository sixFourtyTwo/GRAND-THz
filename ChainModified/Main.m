clear;

p=generate_channel_param_TIV_mp_indoor(0);

channel_type='Gaussian';
corleated_type='uncorleated';
p.SNRdB_list = [0:5:30];
p.trials=100;
tic;

ber_turbo_ZF_hard=zeros(length(p.SNRdB_list),p.trials); %zero forcing - hard
ber_turbo_ZF_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_K_best=zeros(length(p.SNRdB_list),p.trials) %K-best
ber_turbo_K_best_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_ZF_soft=zeros(length(p.SNRdB_list),p.trials); %zero forcing - soft
ber_turbo_ZF_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_ML_hard=zeros(length(p.SNRdB_list),p.trials); %Maximum likelihood - hard
ber_turbo_ML_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_ML_soft=zeros(length(p.SNRdB_list),p.trials); %maximum likelihood - soft
ber_turbo_ML_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_SP_hard=zeros(length(p.SNRdB_list),p.trials); %not sure what SP is, never encoutered
ber_turbo_SP_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_SP_soft=zeros(length(p.SNRdB_list),p.trials);
ber_turbo_SP_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_SSD_hard=zeros(length(p.SNRdB_list),p.trials); %subspace detection - hard
ber_turbo_SSD_hard_Unc=zeros(length(p.SNRdB_list),p.trials);
 
ber_turbo_SSD_soft=zeros(length(p.SNRdB_list),p.trials); %subspace detection - soft
ber_turbo_SSD_soft_Unc=zeros(length(p.SNRdB_list),p.trials);


for t=1:p.trials
    result_BER_ZF_hard = zeros(length(p.SNRdB_list),1); % zero forcing detection (linear)
    result_BER_ZF_hard_Unc = zeros(length(p.SNRdB_list),1);


    result_iter_ZF = zeros(length(p.SNRdB_list),1);
    result_BER_ZF = zeros(length(p.SNRdB_list),1);
    result_BER_ZF_Unc = zeros(length(p.SNRdB_list),1);

    result_BER_K_best= zeros(length(p.SNRdB_list),1); % k_best is a sphere decoding algorithm (detection)
    result_BER_K_best_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_K_best = zeros(length(p.SNRdB_list),1); 
    result_BER_K_best = zeros(length(p.SNRdB_list),1);
    result_BER_K_best_Unc = zeros(length(p.SNRdB_list),1);

    result_BER_ML_hard = zeros(length(p.SNRdB_list),1); % maximum likelihood detection - hard
    result_BER_ML_hard_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_ML = zeros(length(p.SNRdB_list),1); 
    result_BER_ML = zeros(length(p.SNRdB_list),1);
    result_BER_ML_Unc = zeros(length(p.SNRdB_list),1);

    result_BER_SP_hard = zeros(length(p.SNRdB_list),1); % SP thing
    result_BER_SP_hard_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_SP = zeros(length(p.SNRdB_list),1); 
    result_BER_SP = zeros(length(p.SNRdB_list),1);
    result_BER_SP_Unc = zeros(length(p.SNRdB_list),1);

      %new:
    result_BER_SSD_hard = zeros(length(p.SNRdB_list),1); % subspace detection - hard
    result_BER_SSD_hard_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_SSD = zeros(length(p.SNRdB_list),1); 
    result_BER_SSD = zeros(length(p.SNRdB_list),1);
    result_BER_SSD_Unc = zeros(length(p.SNRdB_list),1);

    %% this is where the coding needs to be defined
    % placeholders (Nt = number of tx antennas, Mt number of symbols per
    % antenna)
    svd_turbo = zeros(p.Nt*p.Mt,1);
    svd_turbo2 = zeros(p.Nt*p.Mt,1); 
    svx_turbo = zeros(p.Nt*p.Mt,1);
    svx_turbo2 = zeros(p.Nt*p.Mt,1);
    
    % the channel matrix H
    H = zeros(p.Mr*p.Nr*p.frmLen/p.sq,p.Nt*p.Mt);
    dataBits = step(p.decoder.hCRCGen,randi([0 1],p.dataLen-p.decoder.crclen,1));  %generating random data bits to transmit
    
    intrlvrIndices=randperm(p.encLen); [~, iprm] = sort(intrlvrIndices);  %preparing for interleaving
    encBits_turbo = step(p.encoder, dataBits);  %%% major issue here. The TeraMIMO code only supports turbo coding, not sure what to do
    encBits_turbo2=encBits_turbo;
    encBits_turbo = reshape([encBits_turbo(intrlvrIndices);zeros(p.fillerLen,1)],p.sq,p.frmLen/(p.sq))';
    %%
    
    switch channel_type

        case 'Gaussian'
            H=sqrt(0.5)*(randn(p.Mr*p.Nr,p.Nt*p.Mt)+1i*randn(p.Mr*p.Nr,p.Nt*p.Mt));
            Hmat=repmat(H,p.frmLen/p.sq,1);

         case 'rician'
            K = 3; % Rician K-factor, adjust as needed
            H_LOS = ones(p.Mr * p.Nr, p.Nt * p.Mt); % Line-of-sight (LOS) component
            H_NLOS = sqrt(0.5) * (randn(p.Mr * p.Nr, p.Nt * p.Mt) + 1i * randn(p.Mr * p.Nr, p.Nt * p.Mt)); % NLOS component (Rayleigh fading)
            
            H = sqrt(K / (K + 1)) * H_LOS + sqrt(1 / (K + 1)) * H_NLOS; % Rician channel matrix
            
            Hmat = repmat(H, p.frmLen / p.sq, 1);




        case 'thz'
            [CH_Response, CH_Info] = channel_TIV(p, K_abs);
            Hmat=[];
            if(p.Nsubc>=(p.frmLen/p.sq))
                for i=1:p.frmLen/p.sq
                    Hmat=[Hmat;CH_Response.H(:,:,:,i)];
                end
            else
                last_H=CH_Response.H(:,:,:,p.Nsubc);
                Hmat1=repmat(last_H,(p.frmLen/p.sq)-p.Nsubc,1);
                for i=1:p.Nsubc
                    Hmat=[Hmat;CH_Response.H(:,:,:,i)];
                end
                Hmat=[Hmat;Hmat1];
            end
            
    end

        switch corleated_type
            case 'uncorleated'
                uncorrelated_noise = sqrt(0.5)*(randn(p.Mr*p.Nr,1)+1i*randn(p.Mr*p.Nr,1));
                correlated_noise = uncorrelated_noise;

        end

    y_turbo = zeros(p.Mr*p.Nr*p.frmLen/p.sq,1);

    for k=1:length(p.SNRdB_list)
        N0 = sum(p.Es)*10^(-p.SNRdB_list(k)/10);
        sqrtN0 = sqrt(N0); 
        for i=1:p.frmLen/p.sq   
            ofst=0; 
            for l=1:p.Nt*p.Mt
                svd_turbo(l,1) = bi2de(encBits_turbo(i,(1+ofst):(p.q(l)+ofst)),'left-msb');
                svx_turbo(l,1) = p.syms(l,svd_turbo(l,1)+1); 
                ofst=ofst+p.q(l); 
            end
            y_turbo(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr),1) = Hmat(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr), 1:p.Nt*p.Mt) * svx_turbo + sqrtN0 *correlated_noise;
        end

        %KBEST detection 
        K=2;
        svx_out_K_best=KBEST(p,Hmat,y_turbo,k);
        dataBits_out_K_best_Unc=zeros(784,1);
        for tt=1:length(y_turbo)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_K_best(tt));
            dataBits_out_K_best_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_K_best_Unc2=dataBits_out_K_best_Unc(:,1)-0.5;
        databits_out_K_best=tdec(p.decoder,p.llr_clip*dataBits_out_K_best_Unc2(iprm(1:p.encLen)));
        err_turbo_K_best_Unc=(encBits_turbo2~=dataBits_out_K_best_Unc(iprm(1:p.encLen)));

        err_turbo_K_best = (dataBits~=databits_out_K_best);
        result_BER_K_best_Unc(k) = result_BER_K_best_Unc(k) + sum(err_turbo_K_best_Unc);
        result_BER_K_best(k) = result_BER_K_best(k) + sum(err_turbo_K_best);

        % end of KBEST

        % hard detection 
        
        % ZF hard 
        svx_out_ZF=ZF_hard(p,Hmat,y_turbo);
        dataBits_out_ZF_hard_Unc=zeros(784,1);
        for tt=1:length(y_turbo)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_ZF(tt));
            dataBits_out_ZF_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_ZF_hard_Unc2=dataBits_out_ZF_hard_Unc(:,1)-0.5;
        databits_out_ZF_hard=tdec(p.decoder,p.llr_clip*dataBits_out_ZF_hard_Unc2(iprm(1:p.encLen)));
        err_turbo_ZF_hard_Unc=(encBits_turbo2~=dataBits_out_ZF_hard_Unc(iprm(1:p.encLen)));


        
        temp = dataBits_out_ZF_hard_Unc(iprm(1:p.encLen));
        % disp('Input to Dec:');
        % disp(num2str(temp'));
        % disp('Encoded bits before transmission:');
        % disp(num2str(encBits_turbo2'));
        % disp("data bits:");
        % disp(num2str(dataBits'));
        % disp('Decoded output from Dec:');
        % disp(num2str(databits_out_ZF_hard'));

        err_turbo_ZF_hard = (dataBits~=databits_out_ZF_hard);
        result_BER_ZF_hard_Unc(k) = result_BER_ZF_hard_Unc(k) + sum(err_turbo_ZF_hard_Unc);
        result_BER_ZF_hard(k) = result_BER_ZF_hard(k) + sum(err_turbo_ZF_hard);

        disp(result_BER_ZF_hard(k));
        % ML hard 
        svx_out_ML_temp=ML_hard(p,Hmat,y_turbo); svx_out_ML=svx_out_ML_temp(:);
        dataBits_out_ML_hard_Unc=zeros(784,1);
        for tt=1:length(y_turbo)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_ML(tt));
            dataBits_out_ML_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_ML_hard_Unc2=dataBits_out_ML_hard_Unc(:,1)-0.5;
        databits_out_ML_hard=tdec(p.decoder,p.llr_clip*dataBits_out_ML_hard_Unc2(iprm(1:p.encLen)));
        err_turbo_ML_hard_Unc=(encBits_turbo2~=dataBits_out_ML_hard_Unc(iprm(1:p.encLen)));


        err_turbo_ML_hard = (dataBits~=databits_out_ML_hard);
        result_BER_ML_hard_Unc(k) = result_BER_ML_hard_Unc(k) + sum(err_turbo_ML_hard_Unc);
        result_BER_ML_hard(k) = result_BER_ML_hard(k) + sum(err_turbo_ML_hard);

        %SP hard
        svx_out_SP_temp=sphereDecoderHardType(p,Hmat,y_turbo); svx_out_SP=svx_out_SP_temp(:);
        dataBits_out_SP_hard_Unc=zeros(784,1);
        for tt=1:length(y_turbo)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_SP(tt));
            dataBits_out_SP_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_SP_hard_Unc2=dataBits_out_SP_hard_Unc(:,1)-0.5;
        databits_out_SP_hard=tdec(p.decoder,p.llr_clip*dataBits_out_SP_hard_Unc2(iprm(1:p.encLen)));
        err_turbo_SP_hard_Unc=(encBits_turbo2~=dataBits_out_SP_hard_Unc(iprm(1:p.encLen)));


        err_turbo_SP_hard = (dataBits~=databits_out_SP_hard);
        result_BER_SP_hard_Unc(k) = result_BER_SP_hard_Unc(k) + sum(err_turbo_SP_hard_Unc);
        result_BER_SP_hard(k) = result_BER_SP_hard(k) + sum(err_turbo_SP_hard);

        %SSD HARD:
        svx_out_SSD_temp=SSD_hard(p,Hmat,y_turbo); svx_out_SSD=svx_out_SSD_temp(:);
        dataBits_out_SSD_hard_Unc=zeros(784,1);
        for tt=1:length(y_turbo)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_SSD(tt));
            dataBits_out_SSD_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_SSD_hard_Unc2=dataBits_out_SSD_hard_Unc(:,1)-0.5;
        databits_out_SSD_hard=tdec(p.decoder,p.llr_clip*dataBits_out_SSD_hard_Unc2(iprm(1:p.encLen)));
        err_turbo_SSD_hard_Unc=(encBits_turbo2~=dataBits_out_SSD_hard_Unc(iprm(1:p.encLen)));


        err_turbo_SSD_hard = (dataBits~=databits_out_SSD_hard);
        result_BER_SSD_hard_Unc(k) = result_BER_SSD_hard_Unc(k) + sum(err_turbo_SSD_hard_Unc);
        result_BER_SSD_hard(k) = result_BER_SSD_hard(k) + sum(err_turbo_SSD_hard);
        %end Hard part detection 

    end % SNR loop

        ber_turbo_ZF_soft(:,t) = result_BER_ZF(:,1)/p.dataLen;
        ber_turbo_ZF_soft_Unc(:,t) = result_BER_ZF_Unc(:,1)/p.frmLen;
                
        ber_turbo_ZF_hard_Unc(:,t) = result_BER_ZF_hard_Unc(:,1)/(p.frmLen);
        ber_turbo_ZF_hard(:,t) = result_BER_ZF_hard(:,1)/(p.dataLen);

        ber_turbo_K_best_Unc(:,t) = result_BER_K_best_Unc(:,1)/(p.frmLen);
        ber_turbo_K_best(:,t) = result_BER_K_best(:,1)/(p.dataLen);
        
        ber_turbo_ML_hard_Unc(:,t) = result_BER_ML_hard_Unc(:,1)/(p.frmLen);
        ber_turbo_ML_hard(:,t) = result_BER_ML_hard(:,1)/(p.dataLen);
        
        ber_turbo_ML_soft(:,t) = result_BER_ML(:,1)/p.dataLen;
        ber_turbo_ML_soft_Unc(:,t) = result_BER_ML_Unc(:,1)/p.frmLen;

        ber_turbo_SP_hard_Unc(:,t) = result_BER_SP_hard_Unc(:,1)/(p.frmLen);
        ber_turbo_SP_hard(:,t) = result_BER_SP_hard(:,1)/(p.dataLen);
        
        ber_turbo_SP_soft(:,t) = result_BER_SP(:,1)/p.dataLen;
        ber_turbo_SP_soft_Unc(:,t) = result_BER_SP_Unc(:,1)/p.frmLen;

        % Compute BER for SSD Hard decisions
        ber_turbo_SSD_hard_Unc(:,t) = result_BER_SSD_hard_Unc(:,1)/(p.frmLen); % Uncoded BER for SSD Hard
        ber_turbo_SSD_hard(:,t) = result_BER_SSD_hard(:,1)/(p.dataLen);        % Coded BER for SSD Hard

        % Compute BER for SSD Soft decisions
        %ber_turbo_SSD_soft_Unc(:,t) = result_BER_SSD_soft_Unc(:,1)/(p.frmLen); % Uncoded BER for SSD Soft
        %ber_turbo_SSD_soft(:,t) = result_BER_SSD_soft(:,1)/(p.dataLen);        % Coded BER for SSD Soft

end % trials loop

Sim_Duration=toc;
h = floor(Sim_Duration/3600);
min1 = floor( (Sim_Duration - h*3600) / 60 );
sec = Sim_Duration - h*3600 - min1*60;
fprintf('The simulation time (Modulation Classification) is %d h %d min %f sec',h,min1,sec)
set(0,'defaultlinelinewidth',1.5)
%%
figure(1);

plot6 = semilogy(p.SNRdB_list, sum(ber_turbo_ZF_hard_Unc, 2) / p.trials, 'm--x'); hold on;
% plot10 = semilogy(p.SNRdB_list, sum(ber_turbo_SSD_soft_Unc, 2) / p.trials, 'm--*'); hold on;

% Plot for ZF channel equalization
% plot9 = semilogy(p.SNRdB_list, sum(ber_turbo_K_best_Unc, 2) / p.trials, 'b--x'); hold on;
plot10 = semilogy(p.SNRdB_list, sum(ber_turbo_ZF_hard, 2) / p.trials, 'm--*'); hold on;

plot111 = semilogy(p.SNRdB_list, sum(ber_turbo_K_best_Unc, 2) / p.trials, 'b--x'); hold on;
plot222 = semilogy(p.SNRdB_list, sum(ber_turbo_K_best, 2) / p.trials, 'b--*'); hold on;
plot3 = semilogy(p.SNRdB_list, sum(ber_turbo_ML_hard_Unc, 2) / p.trials, 'b--*'); hold on;
plot4 = semilogy(p.SNRdB_list, sum(ber_turbo_ML_hard, 2) / p.trials, 'b--*'); hold on;
plot5 = semilogy(p.SNRdB_list, sum(ber_turbo_SSD_hard_Unc, 2) / p.trials, 'r--^'); hold on;
plot9 = semilogy(p.SNRdB_list, sum(ber_turbo_SSD_hard, 2) / p.trials, 'r--d'); hold on;
plot13 = semilogy(p.SNRdB_list, sum(ber_turbo_SP_hard_Unc, 2) / p.trials, 'c--^'); hold on;
plot14 = semilogy(p.SNRdB_list, sum(ber_turbo_SP_hard, 2) / p.trials, 'c--d'); hold on;

% Set labels and legends
xlabel('SNR (dB)');
ylabel('BER');
legend([plot6, plot10, plot111, plot222, plot3, plot4, plot13, plot14, plot5, plot9],'ZF Hard Un', 'ZF Hard',    'K-Best Hard Un',    'K-Best Hard ',     'ML Hard unc',   'ML Hard ',    'SP Unc', 'SP', 'ssd hard unc', 'ssd hard');
% legend([plot6,plot9],'SSD soft', 'SSD hard');
grid on;