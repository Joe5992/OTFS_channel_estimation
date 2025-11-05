function [oumiga_tao] = omega_tau(M,l)
%输入子载波数M，坐标数l
%输出omega_tau
oumiga_tao = 0;
 for m=0:M-1
     oumiga_tao=oumiga_tao+exp(1i*2*pi*m*(l)/M);
 end
oumiga_tao=oumiga_tao/M;
end

