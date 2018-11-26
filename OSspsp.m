function gspsp = OSspsp(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table)

%This is a step by step transcription of Table 1 in the Obara-Saika paper
alphaRPQ2 = alpha*RPQ2;
%pi^2.5 = 17.493418327624862;
Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));

g_SSSS_N = Prefactor*vrrBoys_Table_4(alphaRPQ2,Boys_Table); %Orders 0,1,2,3,4

%PSSS integrals of orders 0 to 3
%Order 0
s_px_s_s_0 = RPB(1)*g_SSSS_N(1)+RWP(1)*g_SSSS_N(2);
s_py_s_s_0 = RPB(2)*g_SSSS_N(1)+RWP(2)*g_SSSS_N(2);
s_pz_s_s_0 = RPB(3)*g_SSSS_N(1)+RWP(3)*g_SSSS_N(2);

%Order 1
s_px_s_s_1 = RPB(1)*g_SSSS_N(2)+RWP(1)*g_SSSS_N(3);
s_py_s_s_1 = RPB(2)*g_SSSS_N(2)+RWP(2)*g_SSSS_N(3);
s_pz_s_s_1 = RPB(3)*g_SSSS_N(2)+RWP(3)*g_SSSS_N(3);

%PSPS integrals of order 0 to 2
%Order 0
s_px_s_px_0 = RQD(1)*s_px_s_s_0+RWQ(1)*s_px_s_s_1+1/2/(p+q)*g_SSSS_N(2);
s_px_s_py_0 = RQD(2)*s_px_s_s_0+RWQ(2)*s_px_s_s_1;
s_px_s_pz_0 = RQD(3)*s_px_s_s_0+RWQ(3)*s_px_s_s_1;

s_py_s_px_0 = RQD(1)*s_py_s_s_0+RWQ(1)*s_py_s_s_1;
s_py_s_py_0 = RQD(2)*s_py_s_s_0+RWQ(2)*s_py_s_s_1+1/2/(p+q)*g_SSSS_N(2);
s_py_s_pz_0 = RQD(3)*s_py_s_s_0+RWQ(3)*s_py_s_s_1;

s_pz_s_px_0 = RQD(1)*s_pz_s_s_0+RWQ(1)*s_pz_s_s_1;
s_pz_s_py_0 = RQD(2)*s_pz_s_s_0+RWQ(2)*s_pz_s_s_1;
s_pz_s_pz_0 = RQD(3)*s_pz_s_s_0+RWQ(3)*s_pz_s_s_1+1/2/(p+q)*g_SSSS_N(2);

gspsp = zeros(1,3,1,3);

gspsp(1,1,1,1) = s_px_s_px_0;
gspsp(1,1,1,2) = s_px_s_py_0;
gspsp(1,1,1,3) = s_px_s_pz_0;

gspsp(1,2,1,1) = s_py_s_px_0;
gspsp(1,2,1,2) = s_py_s_py_0;
gspsp(1,2,1,3) = s_py_s_pz_0;

gspsp(1,3,1,1) = s_pz_s_px_0;
gspsp(1,3,1,2) = s_pz_s_py_0;
gspsp(1,3,1,3) = s_pz_s_pz_0;

end