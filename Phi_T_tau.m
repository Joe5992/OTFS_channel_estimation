function Phi_T_tau= Phi_T_tau(M_tau,N_nu,l_max,k_max,x_p,M,N)
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
Phi_T_tau = zeros(M_T*N_T,M_tau*N_nu);
for lie = 0 : M_tau*N_nu - 1
    for k = 0 : N_T - 1
        for l = 0 : M_T - 1
        Phi_T_tau(k*M_T + l + 1,lie + 1) = x_p .* omega_nu(N, -ceil(k_max) + k - gang_k_nu(lie + 1,1)) .* omega_tau_dao(M,l - gang_l_tau(lie + 1,1));
        end
    end
end
end