function [X_after_insert] = Insert_Pilot(X,m_p,n_p,k_max,l_max)
%此处要输入原始的DD域发射矩阵，导频位置，最大时延索引，最大多普勒索引，输出插入导频后的信号矩阵
k_max = ceil(k_max);
for ii = n_p-2*k_max:n_p+2*k_max
    for jj = m_p-l_max:m_p+l_max
        X(ii,jj) = 0;
    end
end
X(n_p,m_p) = 75;%导频的值为10
X_after_insert = X;
end

