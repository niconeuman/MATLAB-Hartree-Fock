function gLSSS_N = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RPC,RPD,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table)

if L1 == 1 %PSSS
    
    g_SSSS_0 = g_SSSS_N(Kab,Kcd,p,q,alpha,RPQ2,0,Boys_Table);
    g_SSSS_1 = g_SSSS_N(Kab,Kcd,p,q,alpha,RPQ2,1,Boys_Table);
    gLSSS_N = RPA*g_SSSS_0-a/p*RPQ*g_SSSS_1;
    
elseif L1 == 2 %DSSS
    g_SSSS_0 = g_SSSS_N(Kab,Kcd,p,q,alpha,RPQ2,0,Boys_Table);
    g_SSSS_1 = g_SSSS_N(Kab,Kcd,p,q,alpha,RPQ2,1,Boys_Table);
    
    gPSSS_N = RPA*g_SSSS_0-a/p*RPQ*g_SSSS_1;
    
end



end