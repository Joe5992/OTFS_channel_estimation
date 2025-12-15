function x=A_v(M_tau,N_nu,l_max,k_max,x_p,M,N,nu_h_t,Sigma_h_t)

p1=Phi_T_nu(M_tau,N_nu,l_max,k_max,x_p,M,N);
x=real(p1'*p1.*((conj(nu_h_t))*(nu_h_t.')+Sigma_h_t.'));
end