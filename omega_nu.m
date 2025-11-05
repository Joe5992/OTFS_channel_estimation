function [oumiga_v] = omega_nu(N,k)
%输出omega_nu，输入时隙数N，坐标数k
%   此处显示详细说明
 oumiga_v= 0;
 for n=0:N-1
    oumiga_v=oumiga_v+exp(-1i*2*pi*n*(k)/N);
 end

oumiga_v=oumiga_v/N;
end
%读者们可以自己算算泰勒展开，跟readme中的（9）式是一样的