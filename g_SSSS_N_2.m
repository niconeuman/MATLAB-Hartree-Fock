function g_SSSS_N = g_SSSS_N_2(Kab,Kcd,p,q,alpha,RPQ2,N,Boys_Table)

%This function calculates the [00|00]^(N) integrals from N = 0 to N = N
%(specified by the calling function)
%This is done to avoid repeating operations.

alphaRPQ2 = alpha*RPQ2;
%pi^2.5 = 17.493418327624862;
Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));
%Check the speed of this code:
% g_SSSS_N = zeros(N+1,1);
% 
% for k = 1:N+1
%     g_SSSS_N(k) = Prefactor*Interpolated_Boys(k-1,alphaRPQ2,Boys_Table);
% end
%Against the speed of this code
g_SSSS_N = Prefactor*Interpolated_Boys_N_3(N,alphaRPQ2,Boys_Table);
end
