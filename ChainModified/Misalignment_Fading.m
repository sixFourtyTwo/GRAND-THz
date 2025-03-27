function [ratio] = Misalignment_Fading(a, w_d, sigma_s, N)
% Misalignment Fading for channel
% ref: 
% [1] Analytical performance assessment of THz wireless systems
% [2] Outage Capacity Optimization for Free-Space Optical Links With Pointing Errors
% hui.chen@kaust.edu.sa

% input: 
% a: radius of the RX effective area;
% w_d: maximum radius of the beam at distance d;
% sigma_s: std of the pointing error displacement at the RX;
% N: number of realizations

% output:
% hp: misalignment fading coefficient

% example: hp = misalignment_fading(0.01, 0.05, 0.01, 5);


% u: intermediate variable
u = sqrt(pi)*a/(sqrt(2)*w_d);
% A_0: fraction of the collected power at r = 0
A_0 = erf(u)^2;
% w_eq2: w_eq^2, related to w_d^2
w_eq2 = w_d^2*sqrt(pi)*erf(u)/(2*u*exp(-u^2));
% gamma: ratio between the equivalent beam width radius and the std of r at RX
gamma = sqrt(w_eq2)/2/sigma_s;
if u > 1
    disp('Beamwidth too Small!');
    hp = zeros(N,1);
elseif gamma > 5        % set a limit for gamma
    hp = ones(N,1)*A_0;
else
    % resolution for generating discrete cdf
    resolution = A_0/1000;
    hp_sample = resolution:resolution:A_0;
    f_hp = gamma^2/A_0^(gamma^2)*hp_sample.^(gamma^2-1); % vector
    % generate CDF
    xCDF=hp_sample;
    yCDF=cumsum(f_hp)/sum(f_hp);
    x1=rand(N,1);
    x1(x1<yCDF(1)) = yCDF(1);  % interpolation resolution
    hp = interp1(yCDF,xCDF,x1);
end
ratio = hp/A_0;
end