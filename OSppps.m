function gppps = OSppps(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table)

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

%Order 2
px_s_s_s_2 = RPA(1)*g_SSSS_N(3)+RWP(1)*g_SSSS_N(4);
py_s_s_s_2 = RPA(2)*g_SSSS_N(3)+RWP(2)*g_SSSS_N(4);
pz_s_s_s_2 = RPA(3)*g_SSSS_N(3)+RWP(3)*g_SSSS_N(4);

%Order 3
px_s_s_s_3 = RPA(1)*g_SSSS_N(4)+RWP(1)*g_SSSS_N(5);
py_s_s_s_3 = RPA(2)*g_SSSS_N(4)+RWP(2)*g_SSSS_N(5);
pz_s_s_s_3 = RPA(3)*g_SSSS_N(4)+RWP(3)*g_SSSS_N(5);

%SPSS integrals of orders 0 to 3
%Order 0
s_px_s_s_0 = RPB(1)*g_SSSS_N(1)+RWP(1)*g_SSSS_N(2);
s_py_s_s_0 = RPB(2)*g_SSSS_N(1)+RWP(2)*g_SSSS_N(2);
s_pz_s_s_0 = RPB(3)*g_SSSS_N(1)+RWP(3)*g_SSSS_N(2);

%Order 1
s_px_s_s_1 = RPB(1)*g_SSSS_N(2)+RWP(1)*g_SSSS_N(3);
s_py_s_s_1 = RPB(2)*g_SSSS_N(2)+RWP(2)*g_SSSS_N(3);
s_pz_s_s_1 = RPB(3)*g_SSSS_N(2)+RWP(3)*g_SSSS_N(3);


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

%Order 1
px_px_s_s_1 = RPB(1)*px_s_s_s_1+RWP(1)*px_s_s_s_2+1/2/p*(g_SSSS_N(2)-alpha/p*g_SSSS_N(3));
px_py_s_s_1 = RPB(2)*px_s_s_s_1+RWP(2)*px_s_s_s_2;
px_pz_s_s_1 = RPB(3)*px_s_s_s_1+RWP(3)*px_s_s_s_2;

py_px_s_s_1 = RPB(1)*py_s_s_s_1+RWP(1)*py_s_s_s_2;
py_py_s_s_1 = RPB(2)*py_s_s_s_1+RWP(2)*py_s_s_s_2+1/2/p*(g_SSSS_N(2)-alpha/p*g_SSSS_N(3));
py_pz_s_s_1 = RPB(3)*py_s_s_s_1+RWP(3)*py_s_s_s_2;

pz_px_s_s_1 = RPB(1)*pz_s_s_s_1+RWP(1)*pz_s_s_s_2;
pz_py_s_s_1 = RPB(2)*pz_s_s_s_1+RWP(2)*pz_s_s_s_2;
pz_pz_s_s_1 = RPB(3)*pz_s_s_s_1+RWP(3)*pz_s_s_s_2+1/2/p*(g_SSSS_N(2)-alpha/p*g_SSSS_N(3));

%PPPS integrals order 0 and 1
%Order 0
px_px_px_s_0 = RQC(1)*px_px_s_s_0+RWQ(1)*px_px_s_s_1+1/2/(p+q)*(s_px_s_s_1+px_s_s_s_1);
px_px_py_s_0 = RQC(2)*px_px_s_s_0+RWQ(2)*px_px_s_s_1;
px_px_pz_s_0 = RQC(3)*px_px_s_s_0+RWQ(3)*px_px_s_s_1;

px_py_px_s_0 = RQC(1)*px_py_s_s_0+RWQ(1)*px_py_s_s_1+1/2/(p+q)*(s_py_s_s_1);
px_py_py_s_0 = RQC(2)*px_py_s_s_0+RWQ(2)*px_py_s_s_1+1/2/(p+q)*(px_s_s_s_1);
px_py_pz_s_0 = RQC(3)*px_py_s_s_0+RWQ(3)*px_py_s_s_1;

px_pz_px_s_0 = RQC(1)*px_pz_s_s_0+RWQ(1)*px_pz_s_s_1+1/2/(p+q)*(s_pz_s_s_1);
px_pz_py_s_0 = RQC(2)*px_pz_s_s_0+RWQ(2)*px_pz_s_s_1;
px_pz_pz_s_0 = RQC(3)*px_pz_s_s_0+RWQ(3)*px_pz_s_s_1+1/2/(p+q)*(px_s_s_s_1);

py_px_px_s_0 = RQC(1)*py_px_s_s_0+RWQ(1)*py_px_s_s_1+1/2/(p+q)*(py_s_s_s_1);
py_px_py_s_0 = RQC(2)*py_px_s_s_0+RWQ(2)*py_px_s_s_1+1/2/(p+q)*(s_px_s_s_1);
py_px_pz_s_0 = RQC(3)*py_px_s_s_0+RWQ(3)*py_px_s_s_1;

py_py_px_s_0 = RQC(1)*py_py_s_s_0+RWQ(1)*py_py_s_s_1;
py_py_py_s_0 = RQC(2)*py_py_s_s_0+RWQ(2)*py_py_s_s_1+1/2/(p+q)*(s_py_s_s_1+py_s_s_s_1);
py_py_pz_s_0 = RQC(3)*py_py_s_s_0+RWQ(3)*py_py_s_s_1;

py_pz_px_s_0 = RQC(1)*py_pz_s_s_0+RWQ(1)*py_pz_s_s_1;
py_pz_py_s_0 = RQC(2)*py_pz_s_s_0+RWQ(2)*py_pz_s_s_1+1/2/(p+q)*(s_pz_s_s_1);
py_pz_pz_s_0 = RQC(3)*py_pz_s_s_0+RWQ(3)*py_pz_s_s_1+1/2/(p+q)*(py_s_s_s_1);

pz_px_px_s_0 = RQC(1)*pz_px_s_s_0+RWQ(1)*pz_px_s_s_1+1/2/(p+q)*(pz_s_s_s_1);
pz_px_py_s_0 = RQC(2)*pz_px_s_s_0+RWQ(2)*pz_px_s_s_1;
pz_px_pz_s_0 = RQC(3)*pz_px_s_s_0+RWQ(3)*pz_px_s_s_1+1/2/(p+q)*(s_px_s_s_1);

pz_py_px_s_0 = RQC(1)*pz_py_s_s_0+RWQ(1)*pz_py_s_s_1;
pz_py_py_s_0 = RQC(2)*pz_py_s_s_0+RWQ(2)*pz_py_s_s_1+1/2/(p+q)*(pz_s_s_s_1);
pz_py_pz_s_0 = RQC(3)*pz_py_s_s_0+RWQ(3)*pz_py_s_s_1+1/2/(p+q)*(s_py_s_s_1);

pz_pz_px_s_0 = RQC(1)*pz_pz_s_s_0+RWQ(1)*pz_pz_s_s_1;
pz_pz_py_s_0 = RQC(2)*pz_pz_s_s_0+RWQ(2)*pz_pz_s_s_1;
pz_pz_pz_s_0 = RQC(3)*pz_pz_s_s_0+RWQ(3)*pz_pz_s_s_1+1/2/(p+q)*(pz_s_s_s_1+s_pz_s_s_1);

gppps = zeros(3,3,3,1);

gppps(1,1,1) = px_px_px_s_0;
gppps(1,1,2) = px_px_py_s_0;
gppps(1,1,3) = px_px_pz_s_0;

gppps(1,2,1) = px_py_px_s_0;
gppps(1,2,2) = px_py_py_s_0;
gppps(1,2,3) = px_py_pz_s_0;

gppps(1,3,1) = px_pz_px_s_0;
gppps(1,3,2) = px_pz_py_s_0;
gppps(1,3,3) = px_pz_pz_s_0;

gppps(2,1,1) = py_px_px_s_0;
gppps(2,1,2) = py_px_py_s_0;
gppps(2,1,3) = py_px_pz_s_0;

gppps(2,2,1) = py_py_px_s_0;
gppps(2,2,2) = py_py_py_s_0;
gppps(2,2,3) = py_py_pz_s_0;

gppps(2,3,1) = py_pz_px_s_0;
gppps(2,3,2) = py_pz_py_s_0;
gppps(2,3,3) = py_pz_pz_s_0;

gppps(3,1,1) = pz_px_px_s_0;
gppps(3,1,2) = pz_px_py_s_0;
gppps(3,1,3) = pz_px_pz_s_0;

gppps(3,2,1) = pz_py_px_s_0;
gppps(3,2,2) = pz_py_py_s_0;
gppps(3,2,3) = pz_py_pz_s_0;

gppps(3,3,1) = pz_pz_px_s_0;
gppps(3,3,2) = pz_pz_py_s_0;
gppps(3,3,3) = pz_pz_pz_s_0;

end