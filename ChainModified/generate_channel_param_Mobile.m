function p  = generate_channel_param_Mobile()
% =========================================================================
% -- Function to generate the required time-invariant/ time-variant mobile
%     channel parameters following the ITU and 3GPP standards
% =========================================================================

% -- Function: p  = generate_channel_param_Mobile()

% -- Input Arguments:
%       NULL (no input is required)

% -- Output Arguments:
%       p: Channel struct that contains the channel parameters

%=================================================

% -- (c) 2021 Simon Tarboush

% -- e-mail: simon.w.tarboush@gmail.com;

% =========================================================================

%% Constants

p.c = 2.9979e8;             % Speed of light in vacuum

%% Define PDP channel type

p.PDP = 'EVA'; % 'CommLett' 'EVA'
% Options: 'SingleTap', 'EVA', 'Veh_A', 'Veh_B', 'EPA', 'Ped_A', 'Ped_B'
%          'ETU', 'TU_A', 'Indoor_HyperLan', 'TU_3GPP', 'RA_3GPP', 'HT_3GPP'

%% Transmission parameters (Fc and BW)

p.Fc = 4e9;                  % Center frequency of the transmission bandwidth (Hz)
p.BW = 15e3*2^3;% 2^9                 % Total channel bandwidth (Hz)
p.Nsubc = p.BW/15e3;                 % Number of sub-bands in each subcarrier

%% MIMO transceiver design

% Number of arrays

p.Mt = 1; % Number of transmitter (row)
p.Nt = 1; % Number of transmitter (column)

p.Mr = 1; % Number of Receiver (row)
p.Nr = 1; % Number of Receiver (column)

p.Qt = p.Mt*p.Nt; % Total number of transmitter
p.Qr = p.Mr*p.Nr; % Total number of receiver

%% Compute Bc, Tc, and TV delay domain response parameters

p.BCApproxMode = 'Most Popular';
% Options are: /'Most Popular'/'Pessimistic'/'Optimistic'/

p.TCApproxMode = 'Geometric Mean';
% Options are: /'Geometric Mean'/'Pessimistic'/'Optimistic'

% % Uncomment for TV channel
p.Vel = 362;      % --> 1.34 KHz vmax               % Velocity of the UE/BS (km/hr)
p.DoppType = 'Int'; % 'Int', 'Frac' (need to add it to MPA)

% p.DopplerSpecShape = 'Jakes';   % Options are: 'Jakes', 'Flat'
% % p.nSamplesperFrame = p.Qt*p.Nsubc*5000;   % Length of input signal in samples for OFDM, DFT-s-OFDM and SC-FDE
% % for FBMC you should change this to
% p.nSamplesperFrame = p.Qt*p.Nsubc*(1+4)*200; % 4 is the overlapping factor
end