function [h, h_DD, Bc, Tau_rms, Tc, fd_max] = generate_mobile_channel(p)

if ~isfield(p,'DopplerSpecShape')
    p.DopplerSpecShape = 'TIV';
    Tc = NaN(1);
    fd_max = NaN(1);
else
    [Tc, fd_max] = get_Tc_FdMax(p.Fc, p.Vel, p.c, p.TCApproxMode);
end

switch p.PDP
    
    case 'SingleTap'
        DEL = 0;
        PWR = 0;
        
        %%%%%% ITU common channel models %%%%%%
    case 'EVA' % Extended Vehicular A model (EVA)
        % Same As MATLAB Multipath Fading Propagation EVA
        % Vehicular A Channel delay-power profile
        % Ref 3GPP TS 36.101 version 14.3.0 Release 14
        % Model: Extended Vehicular A model (EVA)
        % Number of channel taps: 9
        % Delay spread (r.m.s.): 357 ns
        % Maximum excess tap delay (span): 2510 ns
        DEL = 1e-9*[0 30 150 310 370 710 1090 1730 2510]; %relative delay in ns
        PWR = [0 -1.5 -1.4 -3.6 -0.6  -9.1  -7  -12  -16.9];  %relative power in dB
    case 'Veh_A' % Vel 60 km/hr
        DEL = 1e-9*[0 300 700 1100 1700 2500];
        PWR = [0 -1 -9 -10 -15 -20];
    case 'Veh_B' % Vel 60 km/hr
        DEL = 1e-9*[0 300 8900 12900 17100 20000];
        PWR = [-2.5 0 -12.8 -10 -25.2 -16];
    case 'EPA'
        % Extended Pedestrian A model (EPA)
        DEL = 1e-9*[0 30 70 90 110 190 410];
        PWR = [0 -1 -2 -3 -8 -17.2 -20.8];
    case 'Ped_A'
        % Pedestrian A model
        DEL = 1e-9*[0 110 190 410];
        PWR = [0 -9.7 -19.2 -22.8];
    case 'Ped_B'
        % Pedestrian B model
        DEL = 1e-9*[0 200 800 1200 2300 3700];
        PWR = [0 -0.9 -4.9 -8 -7.8 -23.9];
    case 'ETU'
        % Extended Typical Urban model (ETU)
        DEL = 1e-9*[0 50 120 200 230 500 1600 2300 5000];
        PWR = [-1 -1 -1 0 0 0 -3 -5 -7];
    case 'TU_A'
        % Typical Urban A model
        DEL = 1e-6*[0 0.217 0.512 0.514 0.517 0.674 0.882 1.230 1.287 1.311 1.349 1.533 1.535 1.622 1.818 1.836 1.884 1.943 2.048 2.140] ;
        PWR = [-5.7 -7.6 -10.1 -10.2 -10.2 -11.5 -13.4 -16.3 -16.9 -17.1 -17.4 -19.0 -19.0 -19.8 -21.5 -21.6 -22.1 -22.6 -23.5 -24.3]; 
    case 'Indoor_HyperLan'
        % Indoor Channel HyperLan
        DEL = 1e-9*[0 10 20 30 40 50 60 70 80 90 110 140 170 200 240 290 340 390];
        PWR = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18.0 -22.4 -26.7];
        
        %%%%%% 3GPP TR 25.943 channel models %%%%%%
   case 'TU_3GPP'
        % 3GPP Typical Urban
        DEL = 1e-9*[0 217 512 514 517 674 882 1230 1287 1311 1349 1533 1535 1622 1818 1836 1884 1943 2048 2140];
        PWR = [-5.7 -7.6 -10.1 -10.2 -10.2 -11.5 -13.4 -16.3 -16.9 -17.1 -17.4 -19.0 -19.0 -19.8 -21.5 -21.6 -22.1 -22.6 -23.5 -24.3];
    case 'RA_3GPP'
        % 3GPP Rural Area
        DEL = 1e-9*[0 42 101 129 149 245 312 410 469 528];
        PWR = [-5.2 -6.4 -8.4 -9.3 -10 -13.1 -15.3 -18.5 -20.4 -22.4];
    case 'HT_3GPP'
        % 3GPP Hilly Terrain
        DEL = 1e-9*[0 356 441 528 546 609 625 842 916 941 15000 16172 16492 16876 16882 16978 17615 17827 17849 18016];
        PWR = [-3.6 -8.9 -10.2 -11.5 -11.8 -12.7 -13 -16.2 -17.3 -17.7 -17.6 -22.7 -24.1 -25.8 -25.8 -26.2 -29 -29.9 -30 -30.7];
            
    otherwise 
        error('This type isn''t implemented yet !!');
end

Pscal=10.^(PWR/10);% Convert To Scalar Values from dB
% Computes Channel Characteristics
[Bc, Tau_rms] = get_Bc_TauRms(Pscal, DEL, p.BCApproxMode);

Fs = p.BW;
Ts = 1/Fs;

ChD = round(DEL/Ts)+1;
Lh = max(ChD); % channel length in samples
Pg = zeros(1,Lh);
numtaps = length(ChD);

ChP = 10.^(PWR/10);
ChP_norm = ChP/sum(ChP);

for indx_tap = 1:numtaps
    d = ChD(1,indx_tap);
    Pg(1,d)= Pg(1,d)+ChP_norm(1,indx_tap);
end


switch p.DopplerSpecShape
    
    case 'TIV'
        h = zeros(p.Qr, p.Qt, Lh);
        for indx_rx = 1:p.Qr
            for indx_tx = 1:p.Qt
                for indx_ch = 1:Lh
                    h(indx_rx, indx_tx, indx_ch) = sqrt(Pg(indx_ch)/2)*(randn+1i*randn);
                end
            end
        end
        h_DD = [];
    case 'Jakes'
        h = zeros(p.Qr, p.Qt, Lh, p.nSamplesperFrame);
        
        numBLKperTC = p.nSamplesperFrame*Ts*fd_max;
        num_changes = floor(numBLKperTC);%floor(numBLKperTC) to ovoid f=fmax
        
        nu_freq = -num_changes:num_changes;
        nu_len = length(nu_freq);
        h_DD = zeros(p.Qr, p.Qt, Lh, nu_len);
        
        if num_changes == 0
            Sf = 1;% PSD of Doppler Spectrum
        else
            % JAKES shape //when it is CLASSICAL is 1/(1-(f/fd_max)^2)^1/2
            % CLASS from the 3GPP TS 36.101 version 14.3.0 Release 14
            Sf = 1./(pi*fd_max*sqrt(1-(nu_freq/numBLKperTC).^2));
            Sf = Sf/sum(Sf);
        end
        
        for indx_rx = 1:p.Qr
            for indx_tx = 1:p.Qt
                for indx_ch = 1:Lh
                    for indx_nu = 1:nu_len
                        h_DD(indx_rx, indx_tx, indx_ch, indx_nu) = sqrt(Pg(indx_ch)*Sf(indx_nu)/2)*(randn+1i*randn);
                    end
                end
            end
        end
        
        for indx_rx = 1:p.Qr
            for indx_tx = 1:p.Qt
                for indx_ch = 1:Lh
                    h(indx_rx, indx_tx, indx_ch, :) = ...
                        exp(1j*2*pi/p.nSamplesperFrame*(0:p.nSamplesperFrame-1)'*(-num_changes:num_changes))*squeeze(h_DD(indx_rx, indx_tx, indx_ch,:));
                end
            end
        end
    case 'Flat'
        h = zeros(p.Qr, p.Qt, Lh, p.nSamplesperFrame);
        
        numBLKperTC = p.nSamplesperFrame*Ts*fd_max;
        num_changes = floor(numBLKperTC);%floor(numBLKperTC) to ovoid f=fmax
        
        nu_freq = -num_changes:num_changes;
        nu_len = length(nu_freq);
        h_DD = zeros(p.Qr, p.Qt, Lh, nu_len);
        
        if num_changes == 0
            Sf = 1;% PSD of Doppler Spectrum
        else
            Sf = ones(1,nu_len)/(2*fd_max);
            Sf = Sf/sum(Sf);
        end
        
        for indx_rx = 1:p.Qr
            for indx_tx = 1:p.Qt
                for indx_ch = 1:Lh
                    for indx_nu = 1:nu_len
                        h_DD(indx_rx, indx_tx, indx_ch, indx_nu) = sqrt(Pg(indx_ch)*Sf(indx_nu)/2)*(randn+1i*randn);
                    end
                end
            end
        end
        
        for indx_rx = 1:p.Qr
            for indx_tx = 1:p.Qt
                for indx_ch = 1:Lh
                    h(indx_rx, indx_tx, indx_ch, :) = ...
                        exp(1j*2*pi/p.nSamplesperFrame*(0:p.nSamplesperFrame-1)'*(-num_changes:num_changes))*squeeze(h_DD(indx_rx, indx_tx, indx_ch,:));
                end
            end
        end
    otherwise
        error('This type isn''t implemented yet !!');
end

end