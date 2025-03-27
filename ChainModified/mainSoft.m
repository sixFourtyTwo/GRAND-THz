clear;

p=generate_channel_param_TIV_mp_indoor(0);

channel_type='Gaussian';
p.SNRdB_list = [0:5:50];p1.SNRdB_list = [0:5:50];
p.trials=100;
tic;
K_abs = compute_Abs_Coef(p);
ber_turbo_ZF_hard=zeros(length(p.SNRdB_list),p.trials);
ber_turbo_ZF_hard_Unc=zeros(length(p.SNRdB_list),p.trials);

ber_turbo_ZF=zeros(length(p.SNRdB_list),p.trials);
ber_turbo_ZF_Unc=zeros(length(p.SNRdB_list),p.trials);


for t=1:p.trials
    result_BER_ZF_hard = zeros(length(p.SNRdB_list),1);
    result_BER_ZF_hard_Unc = zeros(length(p.SNRdB_list),1);

    result_iter_ZF = zeros(length(p.SNRdB_list),1); 
    result_BER_ZF = zeros(length(p.SNRdB_list),1);
    result_BER_ZF_Unc = zeros(length(p.SNRdB_list),1);


    svd_turbo = zeros(p.Nt*p.Mt,1);svd_turbo2 = zeros(p.Nt*p.Mt,1); svx_turbo = zeros(p.Nt*p.Mt,1);svx_turbo2 = zeros(p.Nt*p.Mt,1);
    H = zeros(p.Mr*p.Nr*p.frmLen/p.sq,p.Nt*p.Mt);
    dataBits = step(p.decoder.hCRCGen,randi([0 1],p.dataLen-p.decoder.crclen,1));  
    
    intrlvrIndices=randperm(p.encLen); [~, iprm] = sort(intrlvrIndices);  
    encBits_turbo = step(p.encoder, dataBits);  
    encBits_turbo2=encBits_turbo;
    encBits_turbo = reshape([encBits_turbo(intrlvrIndices);zeros(p.fillerLen,1)],p.sq,p.frmLen/(p.sq))';
    switch channel_type
        case 'Gaussian'
            H=sqrt(0.5)*(randn(p.Mr*p.Nr,p.Nt*p.Mt)+1i*randn(p.Mr*p.Nr,p.Nt*p.Mt));
            Hmat=repmat(H,p.frmLen/p.sq,1);
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
            n = sqrt(0.5)*(randn(p.Mr*p.Nr,1)+1i*randn(p.Mr*p.Nr,1));
            y_turbo(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr),1) = Hmat(((i-1)*p.Mr*p.Nr+1):(i*p.Mr*p.Nr),1:p.Nt*p.Mt)*svx_turbo+ sqrtN0*n;                      
        end
        

        % hard detection 

        svx_out_ZF=ZF_hard(p,Hmat,y_turbo);
        dataBits_out_ZF_hard_Unc=zeros(784,1);
        for tt=1:length(y_turbo)
            [~,svd_temp] = find(p.syms(1,:)==svx_out_ZF(tt));
            dataBits_out_ZF_hard_Unc((tt-1)*p.q(1)+1:(tt)*p.q(1),1) = de2bi(real(svd_temp-1),p.q(1),'left-msb');
        end
        
        dataBits_out_ZF_hard_Unc2=dataBits_out_ZF_hard_Unc(:,1)-0.5;
        databits_out_ZF_hard=tdec(p.decoder,p.llr_clip*dataBits_out_ZF_hard_Unc2(iprm(1:p.encLen))); %usually where grand goes in hard detection
        err_turbo_ZF_hard_Unc=(encBits_turbo2~=dataBits_out_ZF_hard_Unc(iprm(1:p.encLen)));

        
        err_turbo_ZF_hard = (dataBits~=databits_out_ZF_hard);
        
        result_BER_ZF_hard_Unc(k) = result_BER_ZF_hard_Unc(k) + sum(err_turbo_ZF_hard_Unc);
        
        result_BER_ZF_hard(k) = result_BER_ZF_hard(k) + sum(err_turbo_ZF_hard);



        % soft detection

        LLR1_ext_ZF = ZF_soft(p,Hmat,y_turbo,N0);
        LLR1_ext_ZF = LLR1_ext_ZF';
        LLR2 = LLR1_ext_ZF(:);
        LLR = LLR2(iprm(1:p.encLen)); % remove filler bits; inv-permute %/\%remove inv-permute for BCH

        [databitshat, LLR3, crcpass_ZF, iter] = tdec(p.decoder,LLR); %grand goes here

        result_iter_ZF(k) = result_iter_ZF(k) + iter; 
        tt3=LLR>0;
        err_turbo_ZF_Unc = (encBits_turbo2~=tt3);
        err_turbo_ZF = (dataBits~=databitshat);
        result_BER_ZF(k) = result_BER_ZF(k) + sum(err_turbo_ZF);
        result_BER_ZF_Unc(k) = result_BER_ZF_Unc(k) + sum(err_turbo_ZF_Unc);
        LLRout = [LLR3(1:768,1);zeros(12,1)]; % decoder generates extrinsic llrs
        

    end % SNR loop

    ber_turbo_ZF(:,t) = result_BER_ZF(:,1)/p.dataLen;
    ber_turbo_ZF_Unc(:,t) = result_BER_ZF_Unc(:,1)/p.frmLen;


    ber_turbo_ZF_hard_Unc(:,t) = result_BER_ZF_hard_Unc(:,1)/(p.frmLen);
    
    ber_turbo_ZF_hard(:,t) = result_BER_ZF_hard(:,1)/(p.dataLen);
    
end % trials loop

Sim_Duration=toc;
h = floor(Sim_Duration/3600);
min1 = floor( (Sim_Duration - h*3600) / 60 );
sec = Sim_Duration - h*3600 - min1*60;
fprintf('The simulation time (Modulation Classification) is %d h %d min %f sec',h,min1,sec)
set(0,'defaultlinelinewidth',1.5)

figure(1)
plot1 = semilogy (p.SNRdB_list,sum(ber_turbo_ZF_hard_Unc,2)/p.trials,'k--x');
xlabel('SNR-dB');
ylabel('BER');
hold on
plot2 = semilogy (p.SNRdB_list,sum(ber_turbo_ZF_hard,2)/p.trials,'m--*');
xlabel('SNR-dB');
ylabel('BER');
hold on
plot3 = semilogy (p.SNRdB_list,sum(ber_turbo_ZF_Unc,2)/p.trials,'b--s');
xlabel('SNR-dB');
ylabel('BER');
hold on
plot4 = semilogy (p.SNRdB_list,sum(ber_turbo_ZF,2)/p.trials,'b--*');
xlabel('SNR-dB');
ylabel('BER');
hold on

legend([plot1,plot2,plot3,plot4],'ZF hard Unc','ZF hard','ZF Unc','ZF');
grid on

