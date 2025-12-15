function [oumiga_tao_dao] = omega_tau_dao(M,l)
%输出omega_tau的偏导数，输入M和位置索引l
oumiga_tao_dao = 0;
 for m=0:M-1
     oumiga_tao_dao = oumiga_tao_dao + (1i*2*pi*m / M).*exp(1i*2*pi*m*(l)/M);
 end
oumiga_tao_dao = -oumiga_tao_dao/M;
end

