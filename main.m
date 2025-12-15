%References
%  [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285
%  [R2]. T. Thaj and E. Viterbo, "Low Complexity Iterative Rake Decision Feedback Equalizer for Zero-Padded OTFS Systems," in IEEE Transactions on Vehicular Technology, vol. 69, no. 12, pp. 15606-15622, Dec. 2020, doi: 10.1109/TVT.2020.3044276.
close all
clear all
rng('shuffle')
%% OTFS parameters%%%%%%%%%%
% N: 时隙数量
N = 32;
% M: 子载波的数量
M = 32;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% 导频位置
m_p = M/2; %l_p
n_p = N/2; %k_p
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
P=zeros(N*M,N*M);
for j=1:N
    for i=1:M
    E=zeros(M,N);
    E(i,j)=1;
    P((j-1)*M+1:j*M,(i-1)*N+1:i*N)=E;
    end
end
Im=eye(M);
% OTFS variant: (RZP / RCP / CP / ZP)
variant='CP';

%% delay-Doppler grid symbol placement
if(strcmp(variant,'ZP'))         
    length_ZP = M/16; % ZP length (required only for ZP-OTFS)
    length_CP = 0;
elseif(strcmp(variant,'CP'))
    length_ZP = 0;
    length_CP = 0; % CP length (required only for CP-OTFS) 
else
    length_ZP=0;
    length_CP=0;
end
M_data=M-length_ZP;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;
% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame
SNR_dB = 10:1:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
%% Initializing simulation error count variables
N_fram = 100;
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
current_frame_number=zeros(1,length(SNR_dB));
avg_NMSE = zeros(1,length(SNR_dB)); %SBL的平均误差
avg_NMSE_2 = zeros(1,length(SNR_dB)); %DFT方法的平均误差
avg_NMSE_3 = zeros(1,length(SNR_dB)); %阈值法的平均误差
avg_NMSE_4 = zeros(1,length(SNR_dB)); %线性估计法的平均误差
%% 开始测试
%此时行代表时隙，列代表子载波
for iesn0 = 1:length(SNR_dB)
    NMSE_total = 0 ;
    NMSE_total_2 = 0 ;
    NMSE_total_3 = 0 ;
    NMSE_total_4 = 0 ;
    FACT_frame = 0;
    for ifram = 1:N_fram
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code.        
        % channel model following 3GPP standard
        c=physconst('LightSpeed');
        max_speed=200;  % km/hr
        nu_max = (max_speed*car_fre)/(c);
        delay_taps = [0  1 2 3 4 ];
        Doppler_taps = [-3.25 2.3 0.24 1.14 3.52 ];
        chan_coef = [0.1+0.1i 0.4 0.24+0.15i 0.2 0.15+0.1i];
        L_set=unique(delay_taps);
        l_max=ceil(max(L_set));  
        k_max=nu_max/(1/(N*T));
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        X = Generate_2D_data_grid(N,M,data,data_grid);    
        X = X.'; 
        %插入导频
        X = Insert_Pilot(X,m_p,n_p,k_max,l_max);
        x_p = X(n_p,m_p);
        %% 发送矩阵经过信道%%%%%                       
        Y = DD_Channel_Output(X,N,M,delay_taps,Doppler_taps,chan_coef,sigma_2,iesn0,T);
        %% 下面开始信道估计
        CE=martix_for_CE(Y,l_max,k_max,m_p,n_p);%用于信道估计的矩阵
        YUZHI = 1;
        CEE = CE.';
        %% 稀疏贝叶斯SBL方法
        tilde_h = zeros(length(chan_coef),1);
        M_T = l_max + 1;
        N_T = 2.*ceil(k_max) + 1;
        for ii = 1:length(chan_coef)
            tilde_h(ii,1) = chan_coef(ii).* exp(-1i.*2.*pi.*(Doppler_taps(ii)./(N*T)).*delay_taps(ii));
        end
        Phi_T = zeros(M_T*N_T,length(chan_coef));
        for lie = 0 : length(chan_coef) - 1
            for k = 0 : N_T - 1
                for l = 0 : M_T - 1
                    Phi_T(k*M_T + l + 1,lie + 1) = x_p .* omega_nu(N, k - ceil(k_max) - Doppler_taps(lie + 1)) .* omega_tau(M,l - delay_taps(lie + 1));
                end
            end
        end
        y_T_ideal = Phi_T * tilde_h;
        y_T = reshape(CEE,65,1);%这里如果改动参数需要把65换掉，很抱歉在这图省事写了硬编码
        yibuxilong_tolerance = 0.1;
        max_iteration_T = 20;
        rou = 10^(-1);
        c = 10^(-2);
        d = 10^(-2);
        M_T = l_max + 1;
        N_T = 2.*ceil(k_max) + 1;
        M_tau = 2.*l_max;
        N_nu = 4.*ceil(k_max);
        [hat_h,~,~,iot_tau_t,ka_nu_t] = one_D_off_grid_SBL(yibuxilong_tolerance,max_iteration_T,rou,c,d,y_T,M_T,N_T,M_tau,N_nu,l_max,k_max,x_p,M,N);
        y_T_offSBL = Gang_Phi_T(ka_nu_t,iot_tau_t,M_tau,N_nu,l_max,k_max,x_p,M,N)*(hat_h);
        current_NMSE = (((y_T_ideal-y_T_offSBL)') * (y_T_ideal-y_T_offSBL)) ./ (y_T_ideal' * y_T_ideal);
        disp(['贝叶斯的第',num2str(iesn0),'个信噪比点的第',num2str(ifram),'采样NMSE是',num2str(current_NMSE)]);
        %% 创新点1：基于DFT的低复杂度估计法
        hat_Doppler_DFT = [];%DFT估计法估计多普勒的初始数组
        h_omega_for_chan_coef_est = [];%取p个最大的估计值
        k_ii = [];
        l_ii = [];
        hat_h_aaa = [];
        for ii = min(delay_taps):max(delay_taps)
                p_tau = length(find(delay_taps == ii));
                if p_tau == 0
                else
                    current_path_Doppler_taps = DFT_H_EST(ii,p_tau,CEE);
                    [est_h_omega,Doppler_index_k] = sort(abs(CEE(ii+1,:)));
                    h_omega_for_chan_coef_est = horzcat(h_omega_for_chan_coef_est,CEE(ii+1,Doppler_index_k(end - p_tau + 1:end))./x_p); %选几个最大的h_omega的估计值
                    k_ii = horzcat(k_ii,Doppler_index_k(end - p_tau + 1:end) - (1+ceil(k_max)));
                    l_ii = horzcat(l_ii,ii .* ones(1,p_tau));
                    hat_Doppler_DFT = horzcat(hat_Doppler_DFT,current_path_Doppler_taps);
                end
                
         end
        h_omega_for_chan_coef_est = h_omega_for_chan_coef_est.';
        A = zeros(length(delay_taps));
        for hang = 1:length(delay_taps)
            for lie = 1:length(delay_taps)
                 A(hang,lie) = omega_nu(N,k_ii(hang) - hat_Doppler_DFT(lie)) .* omega_tau(M , l_ii(hang) - l_ii(lie)) .* exp(-1i.*2.*pi.*(hat_Doppler_DFT(lie)./(N*T)).*l_ii(lie));
            end   
        end   
           
        hat_h_dft = (100000000.*A) \ (100000000.*h_omega_for_chan_coef_est);
        hat_h_dft = hat_h_dft.';
        current_NMSE_2 = calculate_NMSE(M,N,l_ii,hat_Doppler_DFT,hat_h_dft,delay_taps,Doppler_taps,chan_coef,T);
        disp(['DFT的第',num2str(iesn0),'个信噪比点的第',num2str(ifram),'采样NMSE是',num2str(current_NMSE_2)])
        %% 阈值估计法
        [Position,Zengyi_in_position] =  Yu_zhi_est(CEE,l_max,k_max,x_p,m_p,M,N,T,YUZHI);
        [~,hat_delay_taps,hat_Doppler_taps,hat_chan_coef] = get_parameters(Position,Zengyi_in_position,1,2,T,N,M);
        current_NMSE_3 = calculate_NMSE(M,N,hat_delay_taps,hat_Doppler_taps,hat_chan_coef,delay_taps,Doppler_taps,chan_coef,T);
        disp(['阈值的第',num2str(iesn0),'个信噪比点的第',num2str(ifram),'采样NMSE是',num2str(current_NMSE_3)])
        %% 线性估计法（受噪声影响极其严重！！！！）
        hat_Doppler_linear = [];%线性估计法估计多普勒的初始数组
        h_omega_for_chan_coef_est = [];%取p个最大的估计值
        k_ii = [];
        l_ii = [];
        hat_h_aaa = [];
        for ii = min(delay_taps):max(delay_taps)
            p_tau = length(find(delay_taps == ii));
            if p_tau == 0
            else
                [current_path_Doppler_taps,~] = linear_est_for_multi_path(CEE,p_tau,ii,k_max,N,x_p,M);
                [est_h_omega,Doppler_index_k] = sort(abs(CEE(ii+1,:)));
                h_omega_for_chan_coef_est = horzcat(h_omega_for_chan_coef_est,CEE(ii+1,Doppler_index_k(end - p_tau + 1:end))./x_p); %选几个最大的h_omega的估计值
                k_ii = horzcat(k_ii,Doppler_index_k(end - p_tau + 1:end) - (1+ceil(k_max)));
                l_ii = horzcat(l_ii,ii .* ones(1,p_tau));
                hat_Doppler_linear = horzcat(hat_Doppler_linear,current_path_Doppler_taps);
            end

        end
        h_omega_for_chan_coef_est = h_omega_for_chan_coef_est.';
        A = zeros(length(delay_taps));
        for hang = 1:length(delay_taps)
            for lie = 1:length(delay_taps)
                A(hang,lie) = omega_nu(N,k_ii(hang) - hat_Doppler_linear(lie)) .* omega_tau(M , l_ii(hang) - l_ii(lie)) .* exp(-1i.*2.*pi.*(hat_Doppler_linear(lie)./(N*T)).*l_ii(lie));
            end
        end

        hat_h_dft = (100000000.*A) \ (100000000.*h_omega_for_chan_coef_est);
        hat_h_dft = hat_h_dft.';
        current_NMSE_4 = calculate_NMSE(M,N,l_ii,hat_Doppler_linear,hat_h_dft,delay_taps,Doppler_taps,chan_coef,T);
        disp(['线性的第',num2str(iesn0),'个信噪比点的第',num2str(ifram),'采样NMSE是',num2str(current_NMSE_4)])
        NMSE_total = NMSE_total + current_NMSE;
        NMSE_total_2 = NMSE_total_2 +  current_NMSE_2;
        NMSE_total_3 = NMSE_total_3 +  current_NMSE_3;
        NMSE_total_4 = NMSE_total_4 +  current_NMSE_4;
    end
    avg_NMSE(iesn0) = NMSE_total ./ N_fram;
    avg_NMSE_2(iesn0) = NMSE_total_2 ./ N_fram;
    avg_NMSE_3(iesn0) = NMSE_total_3 ./ N_fram;
    avg_NMSE_4(iesn0) = NMSE_total_4 ./ N_fram;
    disp(['第',num2str(iesn0),'个信噪比点的平均NMSE是',num2str(avg_NMSE_2(iesn0))])
end
% 自己把这俩平均avg_NMSE数组拿出来画图即可，此处省略画图代码