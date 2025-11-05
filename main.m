%% 设置种子
%dd域信道参考文献：Off-grid channel estimation with sparse bayesian learning for otfs systems
%信道估计算法参考文献：P. Raviteja, K. T. Phan, and Y. Hong, “Embedded pilot-aided channel estimation for otfs in delay–doppler channels,” IEEE Transactions on Vehicular Technology, vol. PP, no. 99, pp. 1–1, 2019.
close all
clear all
rng('shuffle')
c=physconst('LightSpeed');
%% OTFS系统参数%%%%%%%%%%
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
M_data=M;
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
% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 0:2.5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
%% Initializing simulation error count variables
N_fram = 100;
current_frame_number=zeros(1,length(SNR_dB));
avg_NMSE = zeros(1,length(SNR_dB));

%% 开始进攻
%此时行代表时隙，列代表子载波
for iesn0 = 1:length(SNR_dB)
    NMSE_total = 0 ;
    FACT_frame = 0;
    for ifram = 1:N_fram
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code.        
        % channel model following 3GPP standard
        max_speed=200;  % km/hr
        nu_max = (max_speed*car_fre)/(c);
        delay_taps = [0 1 2 3 4 ];%时延索引
        Doppler_taps = [-3.2 3.3 4.1 1.14 1.87 ];%多普勒索引
        chan_coef = [0.1+0.4i 0.4 0.5 1 0.55 ];
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
        %x=reshape(X.',N*M,1);
        %% 发送矩阵经过信道%%%%%                       
        Y = DD_Channel_Output(X,N,M,delay_taps,Doppler_taps,chan_coef,sigma_2,iesn0,T);
        %% 下面开始信道估计
        CE=martix_for_CE(Y,l_max,k_max,m_p,n_p);%用于信道估计的矩阵
        YUZHI = 1;
        CEE = CE.';
        [Position,Zengyi_in_position] =  Yu_zhi_est(CEE,l_max,k_max,x_p,m_p,M,N,T,YUZHI);
        [~,hat_delay_taps,hat_Doppler_taps,hat_chan_coef] = get_parameters(Position,Zengyi_in_position,1,2,T,N,M);
        current_NMSE_1 = calculate_NMSE(M,N,hat_delay_taps,hat_Doppler_taps,hat_chan_coef,delay_taps,Doppler_taps,chan_coef,T);
        NMSE_total = NMSE_total +  current_NMSE_1;
        disp(['第',num2str(iesn0),'个信噪比点的第',num2str(ifram),'采样NMSE是',num2str(current_NMSE_1)])
    end
    avg_NMSE(iesn0) = NMSE_total ./ N_fram;
    disp(['第',num2str(iesn0),'个信噪比点的平均NMSE是',num2str(avg_NMSE(iesn0))])
end
figure(1)
semilogy(SNR_dB,avg_NMSE,'r*-','LineWidth',2,'MarkerSize',8)
legend('Threshold-based channel estimation')
grid on
hold on
xlabel('SNR(dB)')
ylabel('NMSE')