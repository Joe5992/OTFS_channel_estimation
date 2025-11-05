function [hat_h_omega,hat_delay_taps,hat_Doppler_taps,hat_chan_coef] = get_parameters(Position,Zengyi_in_position,k,l,T,N,M)
%输入有路径的位置信息，对应位置的增益，输出估计过后的位置k,l处的h_omega
%   此处显示详细说明
num_of_road = sum(sum(Position));%路径的数目
hat_delay_taps = zeros(1,num_of_road);
hat_Doppler_taps = zeros(1,num_of_road);
hat_chan_coef = zeros(1,num_of_road);
[hangshu,lieshu] = size(Position);
Zengyi_in_position = Zengyi_in_position .* Position;
i = 1;
    for hang = 1:hangshu
        for lie = 1:lieshu
        if Zengyi_in_position(hang,lie) ~= 0
            hat_delay_taps(1,i) = hang - 1;
            hat_Doppler_taps(1,i) = lie - ceil(lieshu/2);
            hat_chan_coef(1,i) = Zengyi_in_position(hang,lie);
            i=i+1;
        end
        end
    end
hat_h_omega = h_omega(k,l,hat_chan_coef,hat_delay_taps,hat_Doppler_taps,T,N,M,num_of_road);
end

