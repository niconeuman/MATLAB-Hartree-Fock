function g_SSSS_N = g_SSSS_N(Kab,Kcd,p,q,alpha,RPQ2,N,Boys_Table)

g_SSSS_N = Kab*Kcd*Interpolated_Boys(N,alpha*RPQ2,Boys_Table)*2*pi^2.5/(p*q*sqrt(p+q));


end