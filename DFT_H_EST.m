function [hat_k] = DFT_H_EST(l_tau,p_tau,CEE)
%DFT_H_EST 宋维健的毕设方案，采用DFT-H校准算法估计多普勒的离网分量
%  输入对应的时延l_tau，时延内的信道径数p_tau，多普勒网格数length(CEE(l_tau,:))，信道估计矩阵CEE
%  输出多普勒索引的精确值 = 多普勒离网分量 + 多普勒在网分量
abs_CEE = abs(CEE(l_tau + 1, :));
[~,index_Doppler] = sort(abs_CEE);%粗略地估计在网分量
hat_k = zeros(1,p_tau);
K_P = length(CEE(l_tau + 1,:));
K_P = ceil(K_P./2);
for ii =1:p_tau
   current_index = index_Doppler(end - ii + 1);
   if current_index == 1 ||current_index == length( CEE(l_tau + 1,:) )
       disp('这次不算数');
       continue
   else
   hat_k(ii) =  (current_index + tan(  pi ./ length( CEE(l_tau + 1,:) ) ) ./ (pi ./ length( CEE(l_tau + 1,:) ) ) .* real( (CEE(l_tau + 1,current_index - 1) -CEE(l_tau + 1,current_index + 1) ) ./ (2.* CEE(l_tau + 1,current_index) - CEE(l_tau + 1,current_index + 1) - CEE(l_tau + 1,current_index - 1))  ) )  - K_P;
   end
end
end

