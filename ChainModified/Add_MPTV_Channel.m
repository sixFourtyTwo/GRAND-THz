function  Y = Add_MPTV_Channel(p, q, h, X)
% p: channel struct
% q: waveform struct

% Add FSC/FF Channel Effects && TV/TIV Effects
% X  => p.NTx x p.Ntot
% h  => p.NRx x p.NTx x L
% p.NTx = p.Mt*p.Nt
% p.NRx = p.Mr*p.Nr
% p.Ntot = p.nSubcarr*(1+p.CP+p.CS)
% L: Channel Length in Delay Domain = Tau RMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NT = p.Qt;
% NR = p.Qr;

% Lx = q.Ntot;

ch_tv = false;
if isfield(p,'DopplerSpecShape')
    ch_tv = true;
end
if isfield(p,'PDP')
    % Mobile channel 3D TIV, 4D TV
    if  ~ch_tv % Time In-Variant Channel
        Y = zeros(p.Qr, q.Ntot);
        for nr = 1:p.Qr
            for nt = 1:p.Qt
                Y(nr,:) = Y(nr,:) + filter(squeeze(h(nr,nt,:)),1,X(nt,:));
            end
        end
    else % Time Variant Channel
        Ch_len = size(h,3);
        Y = zeros(p.Qr, q.Ntot);
        for nr = 1:p.Qr
            for nt = 1:p.Qt
                for indx_samp = 1:Ch_len
                    for indx_cht = 1:Ch_len
                        if (indx_samp-indx_cht+1)>0
                            Y(nr,indx_samp) = Y(nr,indx_samp) + h(nr,nt,indx_cht,indx_samp)*X(nt,indx_samp-indx_cht+1);
                        end
                    end
                end
                for indx_samp = Ch_len+1:q.Ntot
                    for indx_cht = 1:Ch_len
                        Y(nr,indx_samp) = Y(nr,indx_samp) + h(nr,nt,indx_cht,indx_samp)*X(nt,indx_samp-indx_cht+1);
                    end
                end
            end
        end 
    end
else
    % THz channel 4D TIV, 5D TV
    if  ~ch_tv % Time In-Variant Channel
        %======= FSC =======
        Y = zeros(p.Qr, q.Ntot);
        %     S = complex(zeros(1,Lx+L-1));%version 1 CP-OFDM
        % v2 S = complex(zeros(1,q.Ntot));% version 2 CP-OFDM and FBMC-OQAM
        for nr = 1:p.Qr
            % v2 S(:) = 0;
            for nt = 1:p.Qt
                htemp = squeeze(h{nr,nt}).'; % MATRIX VERSION squeeze(h(nr,nt,:)).'
                % S = S + conv(htemp,X(nt,:)); %version 1 CP-OFDM
                % v2 S = S + filter(htemp,1,X(nt,:));
                Y(nr,:) = Y(nr,:) + filter(htemp,1,X(nt,:));
            end
            % (nr,:) = S(1:Lx);%version 1 CP-OFDM
            % v2 Y(nr,:) = S;
        end
    else % Time Variant Channel
        error('Under Development!!! ');
    end
end
end
%======= NO CHANNEL =======
% Y=X;