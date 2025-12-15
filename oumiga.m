function [s,oumiga_v,oumiga_tao]= oumiga(k1,k2,k_vi,l1,l2,l_taoi,N,M)
oumiga_v=0;oumiga_tao=0;

for n=0:N-1
    oumiga_v=oumiga_v+exp(-1i*2*pi*n*(k1-k2-k_vi)/N);
end
oumiga_v=oumiga_v/N;

for m=0:M-1
    oumiga_tao=oumiga_tao+exp(1i*2*pi*m*(l1-l2-l_taoi)/M);
end
oumiga_tao=oumiga_tao/M;


s=oumiga_v*oumiga_tao;
end