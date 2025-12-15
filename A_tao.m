function x=A_tao(M_tau,N_nu,l_max,k_max,x_p,M,N,Sigma_h_t,nu_h_t)

p2=Phi_T_tau(M_tau,N_nu,l_max,k_max,x_p,M,N);
x=real(p2'*p2.*(conj(nu_h_t)*(nu_h_t.')+Sigma_h_t.'));

end