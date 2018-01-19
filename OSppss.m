function gppss = OSppss(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table)

%This is a step by step transcription of Table 1 in the Obara-Saika paper
alphaRPQ2 = alpha*RPQ2;
%pi^2.5 = 17.493418327624862;
Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));

g_SSSS_N = Prefactor*vrrBoys_Table_4(alphaRPQ2,Boys_Table); %Orders 0,1,2,3,4

%PSSS integrals of orders 0 to 3
%Order 0
px_s_s_s_0 = RPA(1)*g_SSSS_N(1)+RWP(1)*g_SSSS_N(2);
py_s_s_s_0 = RPA(2)*g_SSSS_N(1)+RWP(2)*g_SSSS_N(2);
pz_s_s_s_0 = RPA(3)*g_SSSS_N(1)+RWP(3)*g_SSSS_N(2);

%Order 1
px_s_s_s_1 = RPA(1)*g_SSSS_N(2)+RWP(1)*g_SSSS_N(3);
py_s_s_s_1 = RPA(2)*g_SSSS_N(2)+RWP(2)*g_SSSS_N(3);
pz_s_s_s_1 = RPA(3)*g_SSSS_N(2)+RWP(3)*g_SSSS_N(3);

%PPSS integrals of order 0 to 2
%Order 0
px_px_s_s_0 = RPB(1)*px_s_s_s_0+RWP(1)*px_s_s_s_1+1/2/p*(g_SSSS_N(1)-alpha/p*g_SSSS_N(2));
px_py_s_s_0 = RPB(2)*px_s_s_s_0+RWP(2)*px_s_s_s_1;
px_pz_s_s_0 = RPB(3)*px_s_s_s_0+RWP(3)*px_s_s_s_1;

py_px_s_s_0 = RPB(1)*py_s_s_s_0+RWP(1)*py_s_s_s_1;
py_py_s_s_0 = RPB(2)*py_s_s_s_0+RWP(2)*py_s_s_s_1+1/2/p*(g_SSSS_N(1)-alpha/p*g_SSSS_N(2));
py_pz_s_s_0 = RPB(3)*py_s_s_s_0+RWP(3)*py_s_s_s_1;

pz_px_s_s_0 = RPB(1)*pz_s_s_s_0+RWP(1)*pz_s_s_s_1;
pz_py_s_s_0 = RPB(2)*pz_s_s_s_0+RWP(2)*pz_s_s_s_1;
pz_pz_s_s_0 = RPB(3)*pz_s_s_s_0+RWP(3)*pz_s_s_s_1+1/2/p*(g_SSSS_N(1)-alpha/p*g_SSSS_N(2));

gppss = zeros(3,3);
gppss(1,1) = px_px_s_s_0;
gppss(1,2) = px_py_s_s_0;
gppss(1,3) = px_pz_s_s_0;

gppss(2,1) = py_px_s_s_0;
gppss(2,2) = py_py_s_s_0;
gppss(2,3) = py_pz_s_s_0;

gppss(3,1) = pz_px_s_s_0;
gppss(3,2) = pz_py_s_s_0;
gppss(3,3) = pz_pz_s_s_0;

end


