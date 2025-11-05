function [NMSE] = calculate_NMSE(M,N,hat_delay_taps,hat_Doppler_taps,hat_chan_coef,delay_taps,Doppler_taps,chan_coef,T)
%CALCULATE_NMSE 输入M和N，真正的信道信息，估计的信道信息，计算h_omega和hat_h_omega并且输出h_omega的NMSE
%   此处显示详细说明
FenZi = 0;
FenMu = 0;
for k = 1:N
    for l = 1:M
        FenMu = FenMu + (abs(h_omega(k,l,chan_coef,delay_taps,Doppler_taps,T,N,M,length(chan_coef)))).^2;
        FenZi = FenZi + abs((h_omega(k,l,chan_coef,delay_taps,Doppler_taps,T,N,M,length(chan_coef)) - h_omega(k,l,hat_chan_coef,hat_delay_taps,hat_Doppler_taps,T,N,M,length(hat_chan_coef)))).^2;
    end
end
NMSE = FenZi ./ FenMu;
end

