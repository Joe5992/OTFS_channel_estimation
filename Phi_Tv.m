function x= Phi_Tv(XP,M_T,N_T,kp,lp,off_k,off_l,N,M,M_tao,N_v)
x=zeros(M_T*N_T,M_tao*N_v);
for k=0:N_T-1
    for l=0:M_T-1
        for i=0:N_v-1
            for j=0:M_tao-1
                [~,~,c]=oumiga(k,kp,off_k(i*M_tao+j+1),l,lp,off_l(i*M_tao+j+1),N,M);
                [a,~]=bias_oumiga(k,kp,off_k(i*M_tao+j+1),l,lp,off_l(i*M_tao+j+1),N,M);
                x(k*M_T+l+1,i*M_tao+j+1)=XP*c*a;
            end
        end
    end
end
end