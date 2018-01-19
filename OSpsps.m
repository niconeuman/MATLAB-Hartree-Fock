function gpsps = OSpsps(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table)

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

%PSPS integrals of order 0 to 2
%Order 0
px_s_px_s_0 = RQC(1)*px_s_s_s_0+RWQ(1)*px_s_s_s_1+1/2/(p+q)*g_SSSS_N(2);
px_s_py_s_0 = RQC(2)*px_s_s_s_0+RWQ(2)*px_s_s_s_1;
px_s_pz_s_0 = RQC(3)*px_s_s_s_0+RWQ(3)*px_s_s_s_1;

py_s_px_s_0 = RQC(1)*py_s_s_s_0+RWQ(1)*py_s_s_s_1;
py_s_py_s_0 = RQC(2)*py_s_s_s_0+RWQ(2)*py_s_s_s_1+1/2/(p+q)*g_SSSS_N(2);
py_s_pz_s_0 = RQC(3)*py_s_s_s_0+RWQ(3)*py_s_s_s_1;

pz_s_px_s_0 = RQC(1)*pz_s_s_s_0+RWQ(1)*pz_s_s_s_1;
pz_s_py_s_0 = RQC(2)*pz_s_s_s_0+RWQ(2)*pz_s_s_s_1;
pz_s_pz_s_0 = RQC(3)*pz_s_s_s_0+RWQ(3)*pz_s_s_s_1+1/2/(p+q)*g_SSSS_N(2);

gpsps = zeros(3,1,3,1);

gpsps(1,1,1,1) = px_s_px_s_0;
gpsps(1,1,2,1) = px_s_py_s_0;
gpsps(1,1,3,1) = px_s_pz_s_0;

gpsps(2,1,1,1) = py_s_px_s_0;
gpsps(2,1,2,1) = py_s_py_s_0;
gpsps(2,1,3,1) = py_s_pz_s_0;

gpsps(3,1,1,1) = pz_s_px_s_0;
gpsps(3,1,2,1) = pz_s_py_s_0;
gpsps(3,1,3,1) = pz_s_pz_s_0;

end