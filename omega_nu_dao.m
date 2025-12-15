function [oumiga_v_dao] = omega_nu_dao(N,k)
%输出omega_nu的偏导数，输入N和位置索引k
 oumiga_v_dao= 0;
 for n=0:N-1
    oumiga_v_dao=oumiga_v_dao+(2i*pi*n/N)*exp(-1i*2*pi*n*(k)/N);
 end

oumiga_v_dao=oumiga_v_dao/N;
end

