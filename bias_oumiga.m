function [v,tao]=bias_oumiga(k1,k2,k_vi,l1,l2,l_taoi,N,M)
v=0;tao=0;
for n=0:N-1
    v=v+1i*2*pi*n/N*exp(-1i*2*pi*n*(k1-k2-k_vi)/N);
end
v=v/N;
for m=0:M-1
    tao=tao+1i*2*pi*m/M*exp(1i*2*pi*m*(l1-l2-l_taoi)/M);
end
tao=-tao/M;
end