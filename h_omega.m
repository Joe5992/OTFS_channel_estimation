function [OUT_h_omega] = h_omega(k,l,chan_coef,delay_taps,Doppler_taps,T,N,M,P)
%输入参数为h_omega的多普勒坐标和时延坐标、信道每径的增益、对应坐标的omega_nu和omega_tau、时延抽头和多普勒抽头、一个符号的时间T、符号个数N、子载波数M、信道径数
%输出h_omega(k,l)

%   此处显示详细说明
OUT_h_omega = 0;
for i =1:P
OUT_h_omega = OUT_h_omega + chan_coef(i) .* omega_nu(N,k - Doppler_taps(i)) .* omega_tau(M , l - delay_taps(i)) .* exp(-1i.*2.*pi.*(Doppler_taps(i)./(N*T)).*delay_taps(i));
end
end

