function [Gang_Phi_T] = Gang_Phi_T(ka_nu,iot_tau,M_tau,N_nu,l_max,k_max,x_p,M,N)
%输出Gang_Phi_T矩阵，输入多普勒和时延的离网数组ka_nu,iot_tau，多普勒域和延迟域中的虚拟采样网格大小M_tau,N_nu
%输入最大时延和最大多普勒索引l_max,k_max，输入导频位置l_p,k_p，导频本尊x_p，矩阵大小M,N
M_T = l_max + 1;
N_T = 2.*ceil(k_max) + 1;
r_tau = l_max ./ M_tau;%时延分辨率
r_nu = (2.*ceil(k_max)) ./ N_nu;%多普勒分辨率
gang_k_nu = zeros(M_tau*N_nu,1);%初始化在网多普勒
gang_l_tau = zeros(M_tau*N_nu,1);%初始化在网时延
for k_pie_pie = 0:N_nu - 1
    for l_pie_pie = 0:M_tau - 1
        gang_k_nu(k_pie_pie*M_tau + l_pie_pie + 1,1) = k_pie_pie .* r_nu - ceil(k_max);
        gang_l_tau(k_pie_pie*M_tau + l_pie_pie + 1,1) = l_pie_pie .* r_tau;
    end
end
Phi_T = zeros(M_T*N_T,M_tau*N_nu);
Phi_T_nu = zeros(M_T*N_T,M_tau*N_nu);
Phi_T_tau = zeros(M_T*N_T,M_tau*N_nu);

for lie = 0 : M_tau*N_nu - 1
    for k = 0 : N_T - 1
        for l = 0 : M_T - 1
        Phi_T(k*M_T + l + 1,lie + 1) = x_p .* omega_nu(N, -ceil(k_max) + k - gang_k_nu(lie + 1,1)) .* omega_tau(M,l - gang_l_tau(lie + 1,1));
        end
    end
end

for lie = 0 : M_tau*N_nu - 1
    for k = 0 : N_T - 1
        for l = 0 : M_T - 1
        Phi_T_nu(k*M_T + l + 1,lie + 1) = x_p .* omega_nu_dao(N, -ceil(k_max) + k - gang_k_nu(lie + 1,1)) .*omega_tau(M,l - gang_l_tau(lie + 1,1));
        end
    end
end

for lie = 0 : M_tau*N_nu - 1
    for k = 0 : N_T - 1
        for l = 0 : M_T - 1
        Phi_T_tau(k*M_T + l + 1,lie + 1) = x_p .* omega_nu(N, -ceil(k_max) + k - gang_k_nu(lie + 1,1)) .* omega_tau_dao(M,l - gang_l_tau(lie + 1,1));
        end
    end
end
Gang_Phi_T = Phi_T + Phi_T_nu * diag(ka_nu) + Phi_T_tau * diag(iot_tau) ;
end