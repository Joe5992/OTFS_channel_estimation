function [CE_martix] = martix_for_CE(Y_tilda,l_max,k_max,m_p,n_p)
%此处应输入接收到的矩阵，最大时延抽头，最大多普勒抽头，导频的位置，输出用于信道估计的那一部分矩阵
k_max = ceil(k_max);
l_max = ceil(l_max);
CE_martix = zeros(2*k_max+1,l_max+1);
    for k = -k_max:k_max%行索引
        for l = 0:l_max%列索引
            CE_martix(k+1+k_max,l+1) = Y_tilda(n_p+k,m_p+l) ;
        end
    end
end

