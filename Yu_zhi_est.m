function [Position,Zengyi_in_position] = Yu_zhi_est(CEE,l_max,k_max,x_p,l_p,M,N,T,YUZHI)
%YU_ZHI_EST这是阈值估计的函数，输入用于信道估计的矩阵和阈值和导频的值和导频的位置，输出哪个位置可能存在信道，以及对应位置的增益（表示多普勒的K数组、表示时延的L数组）
%   此处显示详细说明
Position = zeros(l_max + 1, 2*ceil(k_max)+1);
Zengyi_in_position = zeros(l_max + 1, 2*ceil(k_max)+1);
z=exp(1i*2*pi/(M*N));
for k = 0:2*ceil(k_max)
    for l = 0:l_max
        if abs(CEE(l+1,k+1)) > YUZHI
            Position(l+1,k+1) = 1;
            Zengyi_in_position(l+1,k+1) = CEE(l+1,k+1)/(x_p.*z^(l_p.*(k+1-(ceil(k_max+1))))).*exp(1i*2*pi*((k+1-(ceil(k_max+1)))./(N*T)).*l);
        end
    end
end
end

