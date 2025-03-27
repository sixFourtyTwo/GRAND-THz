function RxSig = Add_AWGN_Channel(p, TxSig, idx_snr)
% Add AWGN Channel Effects to a Tx Signal

noise = p.sigma(idx_snr)/sqrt(2)*(randn(size(TxSig))+1j*randn(size(TxSig))); % Complex Signal
RxSig = TxSig + noise;

end