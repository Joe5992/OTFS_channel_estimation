function [Y] = DD_Channel_Output(X,N,M,delay_taps,Doppler_taps,chan_coef,sigma_2,iesn0,T)
%输入插入过导频的信道矩阵X,时隙数量N,子载波数量M,信道的真实时延抽头，信道的真实多普勒抽头，每径的信道增益,噪声相关的数据，周期T
% 输出经过信道后的输出矩阵Y，需要写函数来吧h_omega表示出来
Y= zeros(N,M);
for k = 1:N
    for l = 1:M
     %从这里开始是求和
        for k_pie = 1:N
            for l_pie = 1:M
            noise= sqrt(sigma_2(iesn0)/2)*(randn() + 1i*randn());
            Y(k,l)  = Y(k,l) + X(k_pie,l_pie).*h_omega(k-k_pie,l-l_pie,chan_coef,delay_taps,Doppler_taps,T,N,M,length(delay_taps));
            % h_omega由于其有负的索引，故会写成函数
            end
        end
     %求和从这里结束
     %disp(['第',num2str(k),'，',num2str(l),'个元素的值是',num2str(Y(k,l))])
     Y(k,l) = Y(k,l) + noise; %这里在添加噪声
    end
end

end

