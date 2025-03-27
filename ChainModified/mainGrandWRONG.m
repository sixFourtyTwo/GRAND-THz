clear;

p=generate_channel_param_TIV_mp_indoor(0);

channel_type='Gaussian';
corleated_type='uncorleated';
p.SNRdB_list = [0:5:30];
p.trials=100;
tic;

ber_bch_ZF_hard=zeros(length(p.SNRdB_list),p.trials); %zero forcing - hard
ber_bch_ZF_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_K_best=zeros(length(p.SNRdB_list),p.trials) %K-best
ber_bch_K_best_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_ZF_soft=zeros(length(p.SNRdB_list),p.trials); %zero forcing - soft
ber_bch_ZF_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_ML_hard=zeros(length(p.SNRdB_list),p.trials); %Maximum likelihood - hard
ber_bch_ML_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_ML_soft=zeros(length(p.SNRdB_list),p.trials); %maximum likelihood - soft
ber_bch_ML_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_SP_hard=zeros(length(p.SNRdB_list),p.trials); %not sure what SP is, never encoutered
ber_bch_SP_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_SP_soft=zeros(length(p.SNRdB_list),p.trials);
ber_bch_SP_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_SSD_hard=zeros(length(p.SNRdB_list),p.trials); %subspace detection - hard
ber_bch_SSD_hard_Unc=zeros(length(p.SNRdB_list),p.trials);
 
ber_bch_SSD_soft=zeros(length(p.SNRdB_list),p.trials); %subspace detection - soft
ber_bch_SSD_soft_Unc=zeros(length(p.SNRdB_list),p.trials);


parfor t=1:p.trials
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

    % placeholders (Nt = number of tx antennas, Mt number of symbols per
    % antenna)
    svd_bch = zeros(p.Nt*p.Mt,1);
    svd_bch2 = zeros(p.Nt*p.Mt,1);
    svx_bch = zeros(p.Nt*p.Mt,1);
    svx_bch2 = zeros(p.Nt*p.Mt,1);
    
    %% this is where the coding needs to be defined
    % the channel matrix H
    H = zeros(p.Mr*p.Nr*p.frmLen/p.sq,p.Nt*p.Mt);

    bchCode = open('../GRAND-MATLAB-main/CODES/BCH_k_113_n_127.mat'); %obtaining parameters for bch codes
    G = bchCode.G;
    H1 = bchCode.H;
    [np, ~] = size(G);

    dataBits = randi([0 1], np, 1);  %removed CRC, not needed. Removed interleaving

    encBits_bch = mod(dataBits' * G, 2);  %encode using the generator matrix G
    encBits_bch = encBits_bch';

    padding_length = p.encLen - length(encBits_bch);
    encBits_bch_padded = [encBits_bch; zeros(padding_length, 1)];

    encBits_bch2 = encBits_bch_padded;

    encBits_bch_padded = reshape([encBits_bch_padded; zeros(p.fillerLen, 1)], p.sq, p.frmLen / p.sq)';
    
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

    y_bch = zeros(p.Mr*p.Nr*p.frmLen/p.sq,1);

    %this part deals with the transmission, adding noise and so on
    for k=1:length(p.SNRdB_list)
        N0 = sum(p.Es)*10^(-p.SNRdB_list(k)/10);
        sqrtN0 = sqrt(N0);
        for i=1:p.frmLen/p.sq   
            ofst=0; 
            for l=1:p.Nt*p.Mt
                svd_bch(l,1) = bi2de(encBits_bch_padded(i,(1+ofst):(p.q(l)+ofst)),'left-msb');
                svx_bch(l,1) = p.syms(l,svd_bch(l,1)+1); 
                ofst=ofst+p.q(l); 
            end
            y_bch(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr),1) = Hmat(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr), 1:p.Nt*p.Mt) * svx_bch + sqrtN0 *correlated_noise;
        end

        %KBEST detection 
        K=2;
        svx_out_K_best=KBEST(p,Hmat,y_bch,k);
        dataBits_out_K_best_Unc=zeros(784,1);
        for tt=1:length(y_bch)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_K_best(tt));
            dataBits_out_K_best_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb'); %this is of the correct length, but the padding hsa non zero values :(
        end

        dataBits_out_K_best_Unc2=dataBits_out_K_best_Unc(:,1);
        dataBits_out_K_best_Unc2=dataBits_out_K_best_Unc2(1:127);
        [databits_out_K_best_padded,~,~,~] = bin_GRAND(H1,inf,dataBits_out_K_best_Unc2'); %grand goes here
        databits_out_K_best = (databits_out_K_best_padded(1:113))';

        err_bch_K_best_Unc=(encBits_bch2~=dataBits_out_K_best_Unc(1:p.encLen));
        

        err_bch_K_best = (dataBits~=databits_out_K_best);
        result_BER_K_best_Unc(k) = result_BER_K_best_Unc(k) + sum(err_bch_K_best_Unc);
        result_BER_K_best(k) = result_BER_K_best(k) + sum(err_bch_K_best);

        % end of KBEST

        % hard detection 
        
        % ZF hard 
        svx_out_ZF=ZF_hard(p,Hmat,y_bch);
        dataBits_out_ZF_hard_Unc=zeros(784,1);
        for tt=1:length(y_bch)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_ZF(tt));
            dataBits_out_ZF_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end
        
        dataBits_out_ZF_hard_Unc2=dataBits_out_ZF_hard_Unc(:,1);
        dataBits_out_ZF_hard_Unc2 = dataBits_out_ZF_hard_Unc2(1:127);
        [databits_out_ZF_hard_padded,~,nG,~] = bin_GRAND(H1,inf,dataBits_out_ZF_hard_Unc2'); %grand goes here
        databits_out_ZF_hard = databits_out_ZF_hard_padded(1:113)';
        
        err_bch_ZF_hard_Unc=(encBits_bch2~=dataBits_out_ZF_hard_Unc(1:p.encLen));
        
        % %% debug only
        % disp('Input to GRAND:');
        % disp(num2str(dataBits_out_ZF_hard_Unc2'));
        % disp('Encoded bits before transmission:');
        % disp(num2str(encBits_bch2'));
        % disp("DataBits: ");
        % disp(num2str(dataBits'));
        % disp('Decoded output from GRAND:');
        % disp(num2str(databits_out_ZF_hard_padded));
        % 
        %disp("Number of guesses: " + nG);        
        
        
        err_bch_ZF_hard = (dataBits~=databits_out_ZF_hard);
        result_BER_ZF_hard_Unc(k) = result_BER_ZF_hard_Unc(k) + sum(err_bch_ZF_hard_Unc);
        result_BER_ZF_hard(k) = result_BER_ZF_hard(k) + sum(err_bch_ZF_hard);
        
        % ML hard 
        svx_out_ML_temp=ML_hard(p,Hmat,y_bch); svx_out_ML=svx_out_ML_temp(:);
        dataBits_out_ML_hard_Unc=zeros(784,1);
        for tt=1:length(y_bch)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_ML(tt));
            dataBits_out_ML_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end
        
        dataBits_out_ML_hard_Unc2=dataBits_out_ML_hard_Unc(:,1);
        dataBits_out_ML_hard_Unc2=dataBits_out_ML_hard_Unc2(1:127);
        [databits_out_ML_hard_padded,~,~,~] = bin_GRAND(H1,inf,dataBits_out_ML_hard_Unc2'); %grand goes here
        databits_out_ML_hard = databits_out_ML_hard_padded(1:113)';

        err_bch_ML_hard_Unc=(encBits_bch2~=dataBits_out_ML_hard_Unc(1:p.encLen));


        err_bch_ML_hard = (dataBits~=databits_out_ML_hard);
        result_BER_ML_hard_Unc(k) = result_BER_ML_hard_Unc(k) + sum(err_bch_ML_hard_Unc);
        result_BER_ML_hard(k) = result_BER_ML_hard(k) + sum(err_bch_ML_hard);

        %SP hard
        svx_out_SP_temp=sphereDecoderHardType(p,Hmat,y_bch); svx_out_SP=svx_out_SP_temp(:);
        dataBits_out_SP_hard_Unc=zeros(784,1);
        for tt=1:length(y_bch)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_SP(tt));
            dataBits_out_SP_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_SP_hard_Unc2=dataBits_out_SP_hard_Unc(:,1);
        dataBits_out_SP_hard_Unc2=dataBits_out_SP_hard_Unc2(1:127);
        [databits_out_SP_hard_padded,~,~,~] = bin_GRAND(H1,inf,dataBits_out_SP_hard_Unc2'); %grand goes here
        databits_out_SP_hard = databits_out_SP_hard_padded(1:113)';

        err_bch_SP_hard_Unc=(encBits_bch2~=dataBits_out_SP_hard_Unc(1:p.encLen));


        err_bch_SP_hard = (dataBits~=databits_out_SP_hard);
        result_BER_SP_hard_Unc(k) = result_BER_SP_hard_Unc(k) + sum(err_bch_SP_hard_Unc);
        result_BER_SP_hard(k) = result_BER_SP_hard(k) + sum(err_bch_SP_hard);

        %SSD HARD:
        svx_out_SSD_temp=SSD_hard(p,Hmat,y_bch); svx_out_SSD=svx_out_SSD_temp(:);
        dataBits_out_SSD_hard_Unc=zeros(784,1);
        for tt=1:length(y_bch)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_SSD(tt));
            dataBits_out_SSD_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_SSD_hard_Unc2=dataBits_out_SSD_hard_Unc(:,1);
        dataBits_out_SSD_hard_Unc2=dataBits_out_SSD_hard_Unc2(1:127);
        [databits_out_SSD_hard_padded,~,~,~] = bin_GRAND(H1,inf,dataBits_out_SSD_hard_Unc2'); %grand goes here
        databits_out_SSD_hard = databits_out_SSD_hard_padded(1:113)';

        err_bch_SSD_hard_Unc=(encBits_bch2~=dataBits_out_SSD_hard_Unc(1:p.encLen));


        err_bch_SSD_hard = (dataBits~=databits_out_SSD_hard);
        result_BER_SSD_hard_Unc(k) = result_BER_SSD_hard_Unc(k) + sum(err_bch_SSD_hard_Unc);
        result_BER_SSD_hard(k) = result_BER_SSD_hard(k) + sum(err_bch_SSD_hard);
        %end Hard part detection 

    end % SNR loop

        ber_bch_ZF_soft(:,t) = result_BER_ZF(:,1)/p.dataLen;
        ber_bch_ZF_soft_Unc(:,t) = result_BER_ZF_Unc(:,1)/p.frmLen;
                
        ber_bch_ZF_hard_Unc(:,t) = result_BER_ZF_hard_Unc(:,1)/(p.frmLen);
        ber_bch_ZF_hard(:,t) = result_BER_ZF_hard(:,1)/(p.dataLen);

        ber_bch_K_best_Unc(:,t) = result_BER_K_best_Unc(:,1)/(p.frmLen);
        ber_bch_K_best(:,t) = result_BER_K_best(:,1)/(p.dataLen);
        
        ber_bch_ML_hard_Unc(:,t) = result_BER_ML_hard_Unc(:,1)/(p.frmLen);
        ber_bch_ML_hard(:,t) = result_BER_ML_hard(:,1)/(p.dataLen);
        
        ber_bch_ML_soft(:,t) = result_BER_ML(:,1)/p.dataLen;
        ber_bch_ML_soft_Unc(:,t) = result_BER_ML_Unc(:,1)/p.frmLen;

        ber_bch_SP_hard_Unc(:,t) = result_BER_SP_hard_Unc(:,1)/(p.frmLen);
        ber_bch_SP_hard(:,t) = result_BER_SP_hard(:,1)/(p.dataLen);
        
        ber_bch_SP_soft(:,t) = result_BER_SP(:,1)/p.dataLen;
        ber_bch_SP_soft_Unc(:,t) = result_BER_SP_Unc(:,1)/p.frmLen;

        % Compute BER for SSD Hard decisions
        ber_bch_SSD_hard_Unc(:,t) = result_BER_SSD_hard_Unc(:,1)/(p.frmLen); % Uncoded BER for SSD Hard
        ber_bch_SSD_hard(:,t) = result_BER_SSD_hard(:,1)/(p.dataLen);        % Coded BER for SSD Hard

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

plot6 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_hard_Unc, 2) / p.trials, 'm--x'); hold on;
% plot10 = semilogy(p.SNRdB_list, sum(ber_turbo_SSD_soft_Unc, 2) / p.trials, 'm--*'); hold on;

% Plot for ZF channel equalization
% plot9 = semilogy(p.SNRdB_list, sum(ber_turbo_K_best_Unc, 2) / p.trials, 'b--x'); hold on;
plot10 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_hard, 2) / p.trials, 'm-*'); hold on;

plot111 = semilogy(p.SNRdB_list, sum(ber_bch_K_best_Unc, 2) / p.trials, 'g--x'); hold on;
plot222 = semilogy(p.SNRdB_list, sum(ber_bch_K_best, 2) / p.trials, 'g-*'); hold on;
plot3 = semilogy(p.SNRdB_list, sum(ber_bch_ML_hard_Unc, 2) / p.trials, 'b--*'); hold on;
plot4 = semilogy(p.SNRdB_list, sum(ber_bch_ML_hard, 2) / p.trials, 'b-*'); hold on;
plot5 = semilogy(p.SNRdB_list, sum(ber_bch_SSD_hard_Unc, 2) / p.trials, 'r--^'); hold on;
plot9 = semilogy(p.SNRdB_list, sum(ber_bch_SSD_hard, 2) / p.trials, 'r-d'); hold on;
plot13 = semilogy(p.SNRdB_list, sum(ber_bch_SP_hard_Unc, 2) / p.trials, 'c--^'); hold on;
plot14 = semilogy(p.SNRdB_list, sum(ber_bch_SP_hard, 2) / p.trials, 'c-d'); hold on;

% Set labels and legends
xlabel('SNR (dB)');
ylabel('BER');
legend([plot6, plot10, plot111, plot222, plot3, plot4, plot13, plot14, plot5, plot9],'ZF Hard Un', 'ZF Hard',    'K-Best Hard Un',    'K-Best Hard ',     'ML Hard unc',   'ML Hard ',    'SP Unc', 'SP', 'ssd hard unc', 'ssd hard');
% legend([plot6,plot9],'SSD soft', 'SSD hard');
grid on;
