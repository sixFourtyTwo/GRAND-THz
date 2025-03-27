clear;

p=generate_channel_param_TIV_mp_indoor_mod(0);

channel_type='Gaussian';
corleated_type='uncorleated';
p.SNRdB_list = [0:5:30];
p.trials=100;
tic;

ber_bch_ZF_hard=zeros(length(p.SNRdB_list),p.trials); %zero forcing - hard
ber_bch_ZF_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_ZF_soft=zeros(length(p.SNRdB_list),p.trials); %zero forcing - soft
ber_bch_ZF_soft_O1=zeros(length(p.SNRdB_list),p.trials);
ber_bch_ZF_soft_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_ML_hard=zeros(length(p.SNRdB_list),p.trials); %Maximum likelihood - hard
ber_bch_ML_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_bch_ML_soft=zeros(length(p.SNRdB_list),p.trials); %maximum likelihood - soft
ber_bch_ML_soft_O1=zeros(length(p.SNRdB_list),p.trials);
ber_bch_ML_soft_Unc=zeros(length(p.SNRdB_list),p.trials);


parfor t=1:p.trials
    result_BER_ZF_hard = zeros(length(p.SNRdB_list),1); % zero forcing detection (linear)
    result_BER_ZF_hard_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_ZF = zeros(length(p.SNRdB_list),1);

    result_BER_ZF = zeros(length(p.SNRdB_list),1);
    result_BER_ZF_O1 = zeros(length(p.SNRdB_list),1);
    result_BER_ZF_Unc = zeros(length(p.SNRdB_list),1);


    result_BER_ML_hard = zeros(length(p.SNRdB_list),1); % maximum likelihood detection - hard
    result_BER_ML_hard_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_ML = zeros(length(p.SNRdB_list),1);

    result_BER_ML = zeros(length(p.SNRdB_list),1);
    result_BER_ML_O1 = zeros(length(p.SNRdB_list),1);
    result_BER_ML_Unc = zeros(length(p.SNRdB_list),1);


    % placeholders (Nt = number of tx antennas, Mt number of symbols per
    % antenna)
    svd_bch = zeros(p.Nt*p.Mt,1); %%HERE
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

    %encBits_bch_padded = [encBits_bch; zeros(p.fillerLen, 1)];

    encBits_bch2 = encBits_bch;
    
    encBits_bch_padded = reshape([encBits_bch; zeros(p.fillerLen, 1)], p.sq, p.frmLen / p.sq)';
    
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
            if(p.Nsubc>=(p.frmLen/p.sq)) %HERE
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
        N0 = sum(p.Es)*10^(-p.SNRdB_list(k)/10); %HERE
        sqrtN0 = sqrt(N0);
        for i=1:p.frmLen/p.sq   
            ofst=0; 
            for l=1:p.Nt*p.Mt
                svd_bch(l,1) = bi2de(encBits_bch_padded(i,(1+ofst):(p.q(l)+ofst)),'left-msb');
                svx_bch(l,1) = p.syms(l,svd_bch(l,1)+1); %HERE
                ofst=ofst+p.q(l); 
            end
            y_bch(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr),1) = Hmat(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr), 1:p.Nt*p.Mt) * svx_bch + sqrtN0 *correlated_noise;
        end

        % hard detection 
        
        % ZF hard 
        svx_out_ZF=ZF_hard(p,Hmat,y_bch);
        dataBits_out_ZF_hard_Unc=zeros(128,1);
        for tt=1:length(y_bch)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_ZF(tt));
            dataBits_out_ZF_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end

        dataBits_out_ZF_hard_Unc2=dataBits_out_ZF_hard_Unc(1:p.encLen)';
        [databits_out_ZF_hard_padded,~,nG,~] = bin_GRAND(H1,inf,dataBits_out_ZF_hard_Unc2); %grand goes here
        databits_out_ZF_hard = databits_out_ZF_hard_padded(1:p.dataLen)';

        % disp('Input to GRAND:');
        % disp(num2str(dataBits_out_ZF_hard_Unc2));
        % disp('Encoded bits before transmission:');
        % disp(num2str(encBits_bch2'));
        % 
        % disp('Decoded output from GRAND:');
        % disp(num2str(databits_out_ZF_hard_padded'));
        % 
        % disp("Number of guesses: " + nG);

        err_bch_ZF_hard_Unc=(encBits_bch2~=dataBits_out_ZF_hard_Unc(1:p.encLen));
        err_bch_ZF_hard = (dataBits~=databits_out_ZF_hard);
        result_BER_ZF_hard_Unc(k) = result_BER_ZF_hard_Unc(k) + sum(err_bch_ZF_hard_Unc);
        result_BER_ZF_hard(k) = result_BER_ZF_hard(k) + sum(err_bch_ZF_hard);

        %% ZF Soft - ORBGRAND
        LLR1_ext_ZF = ZF_soft(p,Hmat,y_bch,N0);
        LLR1_ext_ZF = LLR1_ext_ZF';
        LLR2 = LLR1_ext_ZF(:);
        LLR = LLR2(1:p.encLen); % remove filler bits; inv-permute %/\%remove inv-permute for BCH
        
        [databitshat,~,nG,~] = bin_ORBGRAND(H1, inf, LLR'); %regular orbgrand
        [databitshat_O1,~,nGO1,~] = bin_ORBGRAND1(H1, inf, LLR'); %orbgrand1
        
        % disp("guesses ORBGRAND: " + nG);
        % disp("guesses ORBGRAND1: " + nGO1);

        databitshat = ~databitshat;
        databitshat_O1 = ~databitshat_O1;

        %result_iter_ZF(k) = result_iter_ZF(k) + iter;
        tt3=LLR>0;
        err_bch_ZF_Unc = (encBits_bch2~=tt3);

        err_bch_ZF = (dataBits~=databitshat(1:p.dataLen)'); %orbgrand
        result_BER_ZF(k) = result_BER_ZF(k) + sum(err_bch_ZF);

        err_bch_ZF_O1 = (dataBits~=databitshat_O1(1:p.dataLen)'); %orbgrand1
        result_BER_ZF_O1(k) = result_BER_ZF_O1(k) + sum(err_bch_ZF_O1);

        result_BER_ZF_Unc(k) = result_BER_ZF_Unc(k) + sum(err_bch_ZF_Unc);
        
        %LLRout = [LLR3(1:768,1);zeros(12,1)]; % decoder generates
        %extrinsic llrs. Assuming this is used for SISO systems
        
        %%

        % % ML hard 
        % svx_out_ML_temp=ML_hard(p,Hmat,y_bch); svx_out_ML=svx_out_ML_temp(:);
        % dataBits_out_ML_hard_Unc=zeros(128,1);
        % for tt=1:length(y_bch)
        %     [~,svd_temp] = find(p.syms(1,:)==svx_out_ML(tt));
        %     dataBits_out_ML_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        % end
        % 
        % dataBits_out_ML_hard_Unc2=dataBits_out_ML_hard_Unc(:,1);
        % [databits_out_ML_hard_padded,~,~,~] = bin_GRAND(H1,inf,dataBits_out_ML_hard_Unc2); %grand goes here
        % databits_out_ML_hard = databits_out_ML_hard_padded(1:113)';
        % 
        % err_bch_ML_hard_Unc=(encBits_bch2~=dataBits_out_ML_hard_Unc);
        % 
        % 
        % err_bch_ML_hard = (dataBits~=databits_out_ML_hard);
        % result_BER_ML_hard_Unc(k) = result_BER_ML_hard_Unc(k) + sum(err_bch_ML_hard_Unc);
        % result_BER_ML_hard(k) = result_BER_ML_hard(k) + sum(err_bch_ML_hard);
        % 
        % %% ML Soft
        % LLR1_ext_ML = ML_soft(p,Hmat,y_bch,N0);
        % LLR1_ext_ML = LLR1_ext_ML';
        % LLR2_ML = LLR1_ext_ML(:);
        % LLR_ML = LLR2_ML(1:127); % remove filler bits; inv-permute %/\%remove inv-permute for BCH
        % 
        % [databitshatpadded_ML,~,~,~] = bin_ORBGRAND(H1, inf, LLR_ML'); %ORBGRAND
        % 
        % [databitshatpadded_ML_O1,~,~,~] = bin_ORBGRAND(H1, inf, LLR_ML'); %ORBGRAND1
        % 
        % databitshat_ML = ~databitshatpadded_ML(1:113)';
        % databitshat_ML_O1 = ~databitshatpadded_ML_O1(1:113)';
        % 
        % %result_iter_ML(k) = result_iter_ML(k) + iter;
        % tt3_ML=LLR2_ML(1:p.encLen)>0;
        % err_bch_ML_Unc = (encBits_bch2~=tt3_ML);
        % 
        % err_bch_ML = (dataBits~=databitshat_ML(1:p.dataLen)); %Orbgrand
        % result_BER_ML(k) = result_BER_ML(k) + sum(err_bch_ML);
        % 
        % err_bch_ML_O1 = (dataBits~=databitshat_ML_O1(1:dataLen)); %Orbgrand1
        % result_BER_ML_O1(k) = result_BER_ML_O1(k) + sum(err_bch_ML_O1);
        % 
        % result_BER_ML_Unc(k) = result_BER_ML_Unc(k) + sum(err_bch_ML_Unc);

    end % SNR loop

        ber_bch_ZF_soft(:,t) = result_BER_ZF(:,1)/p.dataLen;
        ber_bch_ZF_soft_O1(:,t) = result_BER_ZF_O1(:,1)/p.dataLen;
        ber_bch_ZF_soft_Unc(:,t) = result_BER_ZF_Unc(:,1)/p.frmLen;
                
        ber_bch_ZF_hard_Unc(:,t) = result_BER_ZF_hard_Unc(:,1)/(p.frmLen);
        ber_bch_ZF_hard(:,t) = result_BER_ZF_hard(:,1)/(p.dataLen);
        
        ber_bch_ML_hard_Unc(:,t) = result_BER_ML_hard_Unc(:,1)/(p.frmLen);
        ber_bch_ML_hard(:,t) = result_BER_ML_hard(:,1)/(p.dataLen);
        
        ber_bch_ML_soft(:,t) = result_BER_ML(:,1)/p.dataLen;
        ber_bch_ML_soft_O1(:,t) = result_BER_ML_O1(:,1)/p.dataLen;
        ber_bch_ML_soft_Unc(:,t) = result_BER_ML_Unc(:,1)/p.frmLen;

end % trials loop

Sim_Duration=toc;
h = floor(Sim_Duration/3600);
min1 = floor( (Sim_Duration - h*3600) / 60 );
sec = Sim_Duration - h*3600 - min1*60;
fprintf('The simulation time (Modulation Classification) is %d h %d min %f sec',h,min1,sec)
set(0,'defaultlinelinewidth',1.5)
%%
figure(1);

plot11 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_soft, 2) / p.trials, 'r-o'); hold on;
plot12 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_soft_O1, 2) / p.trials, 'g-o'); hold on;
plot13 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_soft_Unc, 2) / p.trials, 'b--o'); hold on;

plot21 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_hard,2) / p.trials, 'r-d'); hold on;
plot22 = semilogy(p.SNRdB_list, sum(ber_bch_ZF_hard_Unc, 2) / p.trials, 'g--d'); hold on;

plot31 = semilogy(p.SNRdB_list, sum(ber_bch_ML_hard, 2) / p.trials, 'r-^'); hold on;
plot32 = semilogy(p.SNRdB_list, sum(ber_bch_ML_hard_Unc, 2) / p.trials, 'g--^'); hold on;

plot41 = semilogy(p.SNRdB_list, sum(ber_bch_ML_soft, 2) / p.trials, 'r-x'); hold on;
plot42 = semilogy(p.SNRdB_list, sum(ber_bch_ML_soft_O1, 2) / p.trials, 'g-x'); hold on;
plot43 = semilogy(p.SNRdB_list, sum(ber_bch_ML_soft_Unc, 2) / p.trials, 'b--x'); hold on;

% Set labels and legends
xlabel('SNR (dB)');
ylabel('BER');
legend([plot11, plot12, plot13, plot21, plot22, plot31, plot32, plot41, plot42, plot43],'ZF-ORBGRAND', 'ZF-ORBGRAND1', 'ZF-ORBGRAND Uncoded', 'ZF-GRAND', 'ZF-GRAND Uncoded', 'ML-GRAND', 'ML-GRAND Uncoded', 'ML-ORBGRAND', 'ML-ORBGRAND1', 'ML-ORBGRAND Uncoded');
% legend([plot6,plot9],'SSD soft', 'SSD hard');
grid on;