function x=b_v(M_tau,N_nu,l_max,k_max,x_p,M,N,nu_h_t,Sigma_h_t,y_T,iot_tau)
x1=Phi_T(M_tau,N_nu,l_max,k_max,x_p,M,N);
x2=Phi_T_tau(M_tau,N_nu,l_max,k_max,x_p,M,N);
x3=Phi_T_nu(M_tau,N_nu,l_max,k_max,x_p,M,N);
x=real(diag(nu_h_t)*(x3.')*conj(y_T)-diag((nu_h_t*(nu_h_t)'+Sigma_h_t)*(x1+x2*diag(iot_tau))' *x3));
end