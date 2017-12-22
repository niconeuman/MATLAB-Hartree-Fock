function g_SSSS_0 = g_SSSS_0(Kab,Kcd,p,q,alpha,RPQ2,Boys_Table)
%This function calculates (SS|SS)^(0) integrals only, in contrast with
%function g_SSSS_N_2, which calculates a series of 0:N order (SS|SS)
%integrals for use in vertical recursion relations.
    %pi^2.5 = 17.493418327624862;
    alphaRPQ2 = alpha*RPQ2;
    Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));
    
    %g_SSSS_0 = Prefactor*Interpolated_Boys(0,alphaRPQ2,Boys_Table);
    tmp = Interpolated_Boys_2(alphaRPQ2,Boys_Table);
    g_SSSS_0 = Prefactor*tmp;

end