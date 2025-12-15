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
c=physconst('LightSpeed');
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
% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 20:2.5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
%% Initializing simulation error count variables
N_fram = 100;
% est_info_bits_MRC=zeros(N_bits_perfram,1);
% err_ber_LMMSE = zeros(1,length(SNR_dB));
% avg_ber_LMMSE = zeros(1,length(SNR_dB));
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
current_frame_number=zeros(1,length(SNR_dB));
avg_NMSE = zeros(1,length(SNR_dB));
avg_NMSE_2 = zeros(1,length(SNR_dB));

%% 开始进攻
%此时行代表时隙，列代表子载波
for iesn0 = 1:length(SNR_dB)
    NMSE_total = 0 ;
    NMSE_total_2 = 0 ;
    FACT_frame = 0;
    for ifram = 1:N_fram
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code.        
        % channel model following 3GPP standard
        max_speed=200;  % km/hr
        nu_max = (max_speed*car_fre)/(c);
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        delay_taps = [0 0 1 2 2 ];
        Doppler_taps = [-3.2 3.3 4.1 1.14 1.87 ];
        chan_coef = [0.1+0.4i 0.4 0.5 1 0.55 ];
        L_set=unique(delay_taps);
        l_max=ceil(max(L_set));  
        k_max=nu_max/(1/(N*T));
%         gs=Gen_time_domain_channel_OTFSvariants(N,M,delay_taps,Doppler_taps,chan_coef,length_CP,variant);
%         G=zeros(N*M,N*M);
%         for n=0:N-1
%             for m=0:M-1
%                 for ell=0:l_max
%                 G(m+n*M+1,n*M+mod(m-ell,M)+1)=gs(ell+1,m+n*M+length_CP+1);
%                 end
%             end
%         end
%         H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');%DD域信道矩阵
%         current_frame_number(iesn0)=ifram;
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
%        noise= sqrt(sigma_2(iesn0)/2)*(randn(size(x)) + 1i*randn(size(x))); 
%         for q=0:N*M-1  
%             n=floor(q/M);
%             m=mod(q,M);
%             for l=0:l_max
%                  r(q+1)=r(q+1)+gs(l+1,m+n*(M+length_CP)+1)*s(n*M+mod(m-l,M)+1);
%             end            
%         end
%         y=H*x;
%         r=P*kron(Im,Fn')*y;
%         r=r+noise;
%         Y_tilda=reshape(r,M,N);     %这一行是经过信道的矩阵
        Y = DD_Channel_Output(X,N,M,delay_taps,Doppler_taps,chan_coef,sigma_2,iesn0,T);
        %% 下面开始信道估计
        CE=martix_for_CE(Y,l_max,k_max,m_p,n_p);%用于信道估计的矩阵
%        %% LMMSE
%         y=reshape(Y.',N*M,1);
%         x_hat1=(H'*H+sigma_2(iesn0).*eye(N*M))^(-1)*(H'*y);
%         x_hat=qamdemod(x_hat1,M_mod,'gray','OutputType','bit');  
%         xx=qamdemod(x,M_mod,'gray','OutputType','bit'); 
%          %% errors count%%%%%
%         errors_LMMSE = sum(xor(x_hat,xx));                
%         err_ber_LMMSE(iesn0) = err_ber_LMMSE(iesn0) + errors_LMMSE; 
        YUZHI = 1;
        CEE = CE.';
        %% 干扰及抗干扰
        random_matrix = zeros(size(CEE));
        random_matrix(5,10) = 10;        
        % 添加干扰信号
        CEE_with_noise = CEE + random_matrix;
       % 已知干扰位置
        interference_pos = [5,10];
        CEE_denoised = cee_recovery_c(CEE_with_noise, interference_pos, n_p, m_p, x_p, M, N);
        delta_angle = angle(CEE(5,10)) - angle(CEE_denoised(5,10));
        %%
        [Position,Zengyi_in_position] =  Yu_zhi_est(CEE,l_max,k_max,x_p,m_p,M,N,T,YUZHI);
        [~,hat_delay_taps,hat_Doppler_taps,hat_chan_coef] = get_parameters(Position,Zengyi_in_position,1,2,T,N,M);
        %[k_d_0,hat_h_0] = linear_est_aft_yuzhi_for_single_path(0,k_max,N,M,x_p,CEE);
        current_NMSE_1 = calculate_NMSE(M,N,hat_delay_taps,hat_Doppler_taps,hat_chan_coef,delay_taps,Doppler_taps,chan_coef,T);
        NMSE_total = NMSE_total +  current_NMSE_1;
        hat_Doppler_DFT = [];%线性估计法估计多普勒的初始数组
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
         if any(abs(hat_Doppler_DFT) > 2.* k_max) || any(abs(hat_h_dft) > 10)
            disp(['DFT的第',num2str(ifram),'次失败了']);
            continue
        else
            FACT_frame = FACT_frame+1;
        end
        hat_h_dft = hat_h_dft.';
        current_NMSE_2 = calculate_NMSE(M,N,l_ii,hat_Doppler_DFT,hat_h_dft,delay_taps,Doppler_taps,chan_coef,T);
        NMSE_total_2 = NMSE_total_2 +  current_NMSE_2;
        disp(['第',num2str(iesn0),'个信噪比点的第',num2str(ifram),'采样NMSE是',num2str(current_NMSE_2)])
    end
    avg_NMSE(iesn0) = NMSE_total ./ N_fram;
    avg_NMSE_2(iesn0) = NMSE_total_2 ./ FACT_frame;
    disp(['第',num2str(iesn0),'个信噪比点的平均NMSE是',num2str(avg_NMSE_2(iesn0))])
end
figure(1)
semilogy(SNR_dB,avg_ber_LMMSE,'r*-','LineWidth',2,'MarkerSize',8)
legend('CP-OTFS')
grid on
hold on
xlabel('SNR(dB)')
ylabel('BER')