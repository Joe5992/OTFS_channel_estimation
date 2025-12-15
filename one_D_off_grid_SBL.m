function [hat_h,hat_l,hat_k,iot_tau_t,ka_nu_t] = one_D_off_grid_SBL(yibuxilong_tolerance,max_iteration_T,rou,c,d,y_T,M_T,N_T,M_tau,N_nu,l_max,k_max,x_p,M,N)
%输出估计增益数组，多普勒数组，时延数组
%输入容忍的容限，最大迭代次数，初始超参数rou，初始超参数c，初始超参数d，导频响应矩阵y_T（你看看reshape一下CEE），输出导频横纵尺寸M_T,N_T，多普勒域和延迟域中的虚拟采样网格大小M_tau,N_nu
%   此处显示详细说明
iteration_count = 1;%迭代计数器
noise_var=y_T'*y_T/100/M_T/N_T;%初始化噪声参数
ka_nu_t = zeros(M_tau*N_nu,1);%多普勒离网数组
iot_tau_t = zeros(M_tau*N_nu,1);%时延离网数组
la_tao_plus = zeros(M_tau*N_nu,1);
ka_v_plus = zeros(M_tau*N_nu,1);
Alpha_t = abs((Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N))' * y_T);%中间那个矩阵Gang_Phi_T(ka_nu,iot_tau)写成了个函数
Beta_0_t = 1 ./ noise_var;
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

while(iteration_count <= max_iteration_T )
    Alpha_t_plus_one = zeros(size(Alpha_t));
    disp(iteration_count);
    Sigma_h_t  = (Beta_0_t .* (Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N))' * (Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N)) + (diag(Alpha_t))^(-1))^(-1); %第t次的方差矩阵(35)
    nu_h_t = Beta_0_t * Sigma_h_t * (Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N))' * y_T;%第t次的均值矩阵(36)
    %超参数更新（Alpha(t + 1)）
    for k_pie_pie = 0:N_nu - 1
        for l_pie_pie = 0:M_tau - 1
            Alpha_t_plus_one(k_pie_pie*M_tau + l_pie_pie + 1) = ( sqrt(1 + 4 .* rou .* ( (abs(nu_h_t(k_pie_pie*M_tau + l_pie_pie + 1))).^2 + Sigma_h_t(k_pie_pie*M_tau + l_pie_pie + 1,k_pie_pie*M_tau + l_pie_pie + 1))) - 1) ./ (2 .* rou);%(44)
        end
    end
    %超参数更新（Beta_0(t + 1)）
    sum=0;
    for k=0:N_nu-1
        for l=0:M_tau-1
            sum=sum+( 1-(Alpha_t(k*M_tau+l+1))^(-1)* Sigma_h_t(k+1,l+1)    );
        end 
    end

    A_Beta_0 = (y_T - Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N) * nu_h_t)'*(y_T - Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N) * nu_h_t) + (Beta_0_t)^(-1) * sum;%(46)
    Beta_0_t_plus_one = (c - 1 + M_T .* N_T)./(d + A_Beta_0);%(45)
    Beta_0_t = Beta_0_t_plus_one;
    %跳出循环(伪代码终止条件)
    if ((Alpha_t_plus_one-Alpha_t)'*(Alpha_t_plus_one-Alpha_t))/(Alpha_t'*Alpha_t) <= yibuxilong_tolerance
        break
    else
        Alpha_t=Alpha_t_plus_one;
    end    
    %更新离网多普勒和时延
     %ka_v(t+1)
    hat_P = ceil((M_T.*N_T) ./ log(M_tau .* N_nu)) - 1;
    [~,index_saved] = sort(Alpha_t);
    A_nu = A_v(M_tau,N_nu,l_max,k_max,x_p,M,N,nu_h_t,Sigma_h_t);
    b_nu = b_v(M_tau,N_nu,l_max,k_max,x_p,M,N,nu_h_t,Sigma_h_t,y_T,iot_tau_t);
    A_v_T = zeros(hat_P,hat_P);
    b_v_T = zeros(hat_P,1);  
    for n = 1:hat_P
        A_v_T(n,n) = A_nu(index_saved(n),index_saved(n));
        b_v_T(n,1) = b_nu(index_saved(n),1);
    end
    ka_nu_T = A_v_T^(-1) * b_v_T;

    for n = 1:hat_P
        ka_v_plus(index_saved(n),1)=ka_nu_T(n,1);
    end

    %%la_tao(t+1)
    b_tau = b_tao(nu_h_t,Sigma_h_t,y_T,M_tau,N_nu,l_max,k_max,x_p,M,N,ka_nu_t);
    A_tau = A_tao(M_tau,N_nu,l_max,k_max,x_p,M,N,Sigma_h_t,nu_h_t);
    A_tao_T = zeros(hat_P,hat_P);
    b_tao_T = zeros(hat_P,1);  
    for n = 1:hat_P
        A_tao_T(n,n) = A_tau(index_saved(n),index_saved(n));
        b_tao_T(n,1) = b_tau(index_saved(n),1);
    end
    iot_tau_T = A_tao_T^(-1) * b_tao_T;
     for n = 1:hat_P
        la_tao_plus(index_saved(n),1)=iot_tau_T(n,1);
    end
   
    iteration_count = iteration_count + 1;
    iot_tau_t=la_tao_plus;
    ka_nu_t=ka_v_plus;
    hat_l = gang_l_tau + iot_tau_t;
    hat_k = gang_k_nu + ka_nu_t;
    hat_h = nu_h_t;
end

end