%参考文献：Efficient Channel Estimation for OTFS Systems in the Presence of Fractional Doppler
%复现场景：一个时延有多个多普勒，利用多项式进行多普勒频移和增益的求解
function [path_Doppler_taps,path_gain_taps] = linear_est_for_multi_path(CEE,p_tau,l,k_max,N,x_p,M)
%LINEAR_EST_FOR_MULTI_PATH 输入导频输出矩阵CEE(h_omega)，已知的径数p_tau，时延抽头l，输出多普勒频移和路径增益
path_Doppler_taps = zeros(1,p_tau);
G_tau = zeros(1,2*ceil(k_max)+1);
for k = -ceil(k_max):ceil(k_max)
    G_tau(k+1+ceil(k_max)) = CEE(l+1,k+1+ceil(k_max))/x_p;
end
[~,Doppler_index_k] = sort(abs(G_tau));%选最大的多普勒抽头
G_i = zeros(1,2*p_tau);
for ii = 0:2*p_tau -1
    G_i(ii + 1) = G_tau(Doppler_index_k(end - ii));%给最大的几个赋值
end
W = exp(-1i.*2.*pi./N);
bf_W = zeros(2*p_tau,2*p_tau);
for hang = 1:2*p_tau
    for lie = 1:p_tau
       bf_W(hang,lie) = W^((Doppler_index_k(end - hang + 1) - 1 - ceil(k_max)).* (lie - 1));
    end
end
for hang = 1:2*p_tau
    for lie = p_tau + 1:2*p_tau
       bf_W(hang,lie) = -G_i(hang) .* W^((Doppler_index_k(end - hang + 1) - 1 - ceil(k_max)).* (lie - 1 - p_tau)) .* W.^(Doppler_index_k(end - hang + 1) - 1 - ceil(k_max));
    end
end
Theta = (100000000.*bf_W) \ (100000000.*G_i.') ;
bf_b = zeros(1 , p_tau + 1);
bf_b(1) = 1;
bf_b(2 : p_tau + 1) = Theta(p_tau + 1:2*p_tau);
bf_z = roots(bf_b);%（26）的解集
for ii = 1 : p_tau 
    path_Doppler_taps(ii) = angle(bf_z(ii)) .* N ./ (2 .* pi);
end
%下面进行增益估计
bf_A = zeros(p_tau,p_tau);
bf_G = zeros(p_tau,1);
for hang = 1:p_tau
    for lie = 1:p_tau
       k_hang = Doppler_index_k(end - hang + 1) - 1 - ceil(k_max);
       z_lie = bf_z(lie);
       bf_A(hang,lie) = ((1 - z_lie.^N).* z_lie.^(-l ./ M)) / (N .* (1 - W.^(k_hang) .* z_lie));
       bf_G(hang,1) = G_i(hang);
    end
end
path_gain_taps = (eye(p_tau)/bf_A )* bf_G;
path_gain_taps = path_gain_taps.';
end

