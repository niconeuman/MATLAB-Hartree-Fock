function gdsds = OSdsds(RPAValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues)

%I am changing these functions so that all the primitive shell specific operations are done in a loop in shellOS
%Then the OSabcd functions work in a vectorized manner

s_s_s_s_0 = gSSSSNValues(:,1);
s_s_s_s_1 = gSSSSNValues(:,2);
s_s_s_s_2 = gSSSSNValues(:,3);
s_s_s_s_3 = gSSSSNValues(:,4);
s_s_s_s_4 = gSSSSNValues(:,5);

PAx = RPAValues(:,1);
PAy = RPAValues(:,2);
PAz = RPAValues(:,3);

QCx = RQCValues(:,1);
QCy = RQCValues(:,2);
QCz = RQCValues(:,3);

WPx = RWPValues(:,1);
WPy = RWPValues(:,2);
WPz = RWPValues(:,3);

WQx = RWQValues(:,1);
WQy = RWQValues(:,2);
WQz = RWQValues(:,3);

oo2p = 0.5./pValues;
oo2q = 0.5./qValues;
qoppq = qValues./ppqValues;
oo2pq = 0.5./ppqValues;

px_s_s_s_0 = PAx*s_s_s_s_0 + WPx*s_s_s_s_1
py_s_s_s_0 = PAy*s_s_s_s_0 + WPy*s_s_s_s_1
pz_s_s_s_0 = PAz*s_s_s_s_0 + WPz*s_s_s_s_1

px_s_s_s_1 = PAx*s_s_s_s_1 + WPx*s_s_s_s_2
py_s_s_s_1 = PAy*s_s_s_s_1 + WPy*s_s_s_s_2
pz_s_s_s_1 = PAz*s_s_s_s_1 + WPz*s_s_s_s_2

px_s_s_s_2 = PAx*s_s_s_s_2 + WPx*s_s_s_s_3
py_s_s_s_2 = PAy*s_s_s_s_2 + WPy*s_s_s_s_3
pz_s_s_s_2 = PAz*s_s_s_s_2 + WPz*s_s_s_s_3

px_s_s_s_3 = PAx*s_s_s_s_3 + WPx*s_s_s_s_4
py_s_s_s_3 = PAy*s_s_s_s_3 + WPy*s_s_s_s_4
pz_s_s_s_3 = PAz*s_s_s_s_3 + WPz*s_s_s_s_4

px_s_px_s_0 = QCx*px_s_s_s_0 + WQx*px_s_s_s_1
px_s_py_s_0 = QCy*px_s_s_s_0 + WQy*px_s_s_s_1
px_s_pz_s_0 = QCz*px_s_s_s_0 + WQz*px_s_s_s_1
py_s_px_s_0 = QCx*py_s_s_s_0 + WQx*py_s_s_s_1
py_s_py_s_0 = QCy*py_s_s_s_0 + WQy*py_s_s_s_1
py_s_pz_s_0 = QCz*py_s_s_s_0 + WQz*py_s_s_s_1
pz_s_px_s_0 = QCx*pz_s_s_s_0 + WQx*pz_s_s_s_1
pz_s_py_s_0 = QCy*pz_s_s_s_0 + WQy*pz_s_s_s_1
pz_s_pz_s_0 = QCz*pz_s_s_s_0 + WQz*pz_s_s_s_1

px_s_px_s_1 = QCx*px_s_s_s_1 + WQx*px_s_s_s_2
px_s_py_s_1 = QCy*px_s_s_s_1 + WQy*px_s_s_s_2
px_s_pz_s_1 = QCz*px_s_s_s_1 + WQz*px_s_s_s_2
py_s_px_s_1 = QCx*py_s_s_s_1 + WQx*py_s_s_s_2
py_s_py_s_1 = QCy*py_s_s_s_1 + WQy*py_s_s_s_2
py_s_pz_s_1 = QCz*py_s_s_s_1 + WQz*py_s_s_s_2
pz_s_px_s_1 = QCx*pz_s_s_s_1 + WQx*pz_s_s_s_2
pz_s_py_s_1 = QCy*pz_s_s_s_1 + WQy*pz_s_s_s_2
pz_s_pz_s_1 = QCz*pz_s_s_s_1 + WQz*pz_s_s_s_2

px_s_px_s_2 = QCx*px_s_s_s_2 + WQx*px_s_s_s_3
px_s_py_s_2 = QCy*px_s_s_s_2 + WQy*px_s_s_s_3
px_s_pz_s_2 = QCz*px_s_s_s_2 + WQz*px_s_s_s_3
py_s_px_s_2 = QCx*py_s_s_s_2 + WQx*py_s_s_s_3
py_s_py_s_2 = QCy*py_s_s_s_2 + WQy*py_s_s_s_3
py_s_pz_s_2 = QCz*py_s_s_s_2 + WQz*py_s_s_s_3
pz_s_px_s_2 = QCx*pz_s_s_s_2 + WQx*pz_s_s_s_3
pz_s_py_s_2 = QCy*pz_s_s_s_2 + WQy*pz_s_s_s_3
pz_s_pz_s_2 = QCz*pz_s_s_s_2 + WQz*pz_s_s_s_3

px_s_dxx_s_0 = QCx*px_s_px_s_0 + WQx*px_s_px_s_1 + 1*oo2q*(px_s_s_s_0 - qoppq*px_s_s_s_1)
px_s_dxy_s_0 = QCx*px_s_py_s_0 + WQx*px_s_py_s_1
px_s_dxz_s_0 = QCx*px_s_pz_s_0 + WQx*px_s_pz_s_1
px_s_dyy_s_0 = QCy*px_s_py_s_0 + WQy*px_s_py_s_1 + 1*oo2q*(px_s_s_s_0 - qoppq*px_s_s_s_1)
px_s_dyz_s_0 = QCy*px_s_pz_s_0 + WQy*px_s_pz_s_1
px_s_dzz_s_0 = QCz*px_s_pz_s_0 + WQz*px_s_pz_s_1 + 1*oo2q*(px_s_s_s_0 - qoppq*px_s_s_s_1)
py_s_dxx_s_0 = QCx*py_s_px_s_0 + WQx*py_s_px_s_1 + 1*oo2q*(py_s_s_s_0 - qoppq*py_s_s_s_1)
py_s_dxy_s_0 = QCx*py_s_py_s_0 + WQx*py_s_py_s_1
py_s_dxz_s_0 = QCx*py_s_pz_s_0 + WQx*py_s_pz_s_1
py_s_dyy_s_0 = QCy*py_s_py_s_0 + WQy*py_s_py_s_1 + 1*oo2q*(py_s_s_s_0 - qoppq*py_s_s_s_1)
py_s_dyz_s_0 = QCy*py_s_pz_s_0 + WQy*py_s_pz_s_1
py_s_dzz_s_0 = QCz*py_s_pz_s_0 + WQz*py_s_pz_s_1 + 1*oo2q*(py_s_s_s_0 - qoppq*py_s_s_s_1)
pz_s_dxx_s_0 = QCx*pz_s_px_s_0 + WQx*pz_s_px_s_1 + 1*oo2q*(pz_s_s_s_0 - qoppq*pz_s_s_s_1)
pz_s_dxy_s_0 = QCx*pz_s_py_s_0 + WQx*pz_s_py_s_1
pz_s_dxz_s_0 = QCx*pz_s_pz_s_0 + WQx*pz_s_pz_s_1
pz_s_dyy_s_0 = QCy*pz_s_py_s_0 + WQy*pz_s_py_s_1 + 1*oo2q*(pz_s_s_s_0 - qoppq*pz_s_s_s_1)
pz_s_dyz_s_0 = QCy*pz_s_pz_s_0 + WQy*pz_s_pz_s_1
pz_s_dzz_s_0 = QCz*pz_s_pz_s_0 + WQz*pz_s_pz_s_1 + 1*oo2q*(pz_s_s_s_0 - qoppq*pz_s_s_s_1)

px_s_dxx_s_1 = QCx*px_s_px_s_1 + WQx*px_s_px_s_2 + 1*oo2q*(px_s_s_s_1 - qoppq*px_s_s_s_2)
px_s_dxy_s_1 = QCx*px_s_py_s_1 + WQx*px_s_py_s_2
px_s_dxz_s_1 = QCx*px_s_pz_s_1 + WQx*px_s_pz_s_2
px_s_dyy_s_1 = QCy*px_s_py_s_1 + WQy*px_s_py_s_2 + 1*oo2q*(px_s_s_s_1 - qoppq*px_s_s_s_2)
px_s_dyz_s_1 = QCy*px_s_pz_s_1 + WQy*px_s_pz_s_2
px_s_dzz_s_1 = QCz*px_s_pz_s_1 + WQz*px_s_pz_s_2 + 1*oo2q*(px_s_s_s_1 - qoppq*px_s_s_s_2)
py_s_dxx_s_1 = QCx*py_s_px_s_1 + WQx*py_s_px_s_2 + 1*oo2q*(py_s_s_s_1 - qoppq*py_s_s_s_2)
py_s_dxy_s_1 = QCx*py_s_py_s_1 + WQx*py_s_py_s_2
py_s_dxz_s_1 = QCx*py_s_pz_s_1 + WQx*py_s_pz_s_2
py_s_dyy_s_1 = QCy*py_s_py_s_1 + WQy*py_s_py_s_2 + 1*oo2q*(py_s_s_s_1 - qoppq*py_s_s_s_2)
py_s_dyz_s_1 = QCy*py_s_pz_s_1 + WQy*py_s_pz_s_2
py_s_dzz_s_1 = QCz*py_s_pz_s_1 + WQz*py_s_pz_s_2 + 1*oo2q*(py_s_s_s_1 - qoppq*py_s_s_s_2)
pz_s_dxx_s_1 = QCx*pz_s_px_s_1 + WQx*pz_s_px_s_2 + 1*oo2q*(pz_s_s_s_1 - qoppq*pz_s_s_s_2)
pz_s_dxy_s_1 = QCx*pz_s_py_s_1 + WQx*pz_s_py_s_2
pz_s_dxz_s_1 = QCx*pz_s_pz_s_1 + WQx*pz_s_pz_s_2
pz_s_dyy_s_1 = QCy*pz_s_py_s_1 + WQy*pz_s_py_s_2 + 1*oo2q*(pz_s_s_s_1 - qoppq*pz_s_s_s_2)
pz_s_dyz_s_1 = QCy*pz_s_pz_s_1 + WQy*pz_s_pz_s_2
pz_s_dzz_s_1 = QCz*pz_s_pz_s_1 + WQz*pz_s_pz_s_2 + 1*oo2q*(pz_s_s_s_1 - qoppq*pz_s_s_s_2)

dxx_s_s_s_0 = PAx*px_s_s_s_0 + WPx*px_s_s_s_1 + 1*oo2p*(s_s_s_s_0 - qoppq*s_s_s_s_1)
dxy_s_s_s_0 = PAx*py_s_s_s_0 + WPx*py_s_s_s_1
dxz_s_s_s_0 = PAx*pz_s_s_s_0 + WPx*pz_s_s_s_1
dyy_s_s_s_0 = PAy*py_s_s_s_0 + WPy*py_s_s_s_1 + 1*oo2p*(s_s_s_s_0 - qoppq*s_s_s_s_1)
dyz_s_s_s_0 = PAy*pz_s_s_s_0 + WPy*pz_s_s_s_1
dzz_s_s_s_0 = PAz*pz_s_s_s_0 + WPz*pz_s_s_s_1 + 1*oo2p*(s_s_s_s_0 - qoppq*s_s_s_s_1)

dxx_s_s_s_1 = PAx*px_s_s_s_1 + WPx*px_s_s_s_2 + 1*oo2p*(s_s_s_s_1 - qoppq*s_s_s_s_2)
dxy_s_s_s_1 = PAx*py_s_s_s_1 + WPx*py_s_s_s_2
dxz_s_s_s_1 = PAx*pz_s_s_s_1 + WPx*pz_s_s_s_2
dyy_s_s_s_1 = PAy*py_s_s_s_1 + WPy*py_s_s_s_2 + 1*oo2p*(s_s_s_s_1 - qoppq*s_s_s_s_2)
dyz_s_s_s_1 = PAy*pz_s_s_s_1 + WPy*pz_s_s_s_2
dzz_s_s_s_1 = PAz*pz_s_s_s_1 + WPz*pz_s_s_s_2 + 1*oo2p*(s_s_s_s_1 - qoppq*s_s_s_s_2)

dxx_s_s_s_2 = PAx*px_s_s_s_2 + WPx*px_s_s_s_3 + 1*oo2p*(s_s_s_s_2 - qoppq*s_s_s_s_3)
dxy_s_s_s_2 = PAx*py_s_s_s_2 + WPx*py_s_s_s_3
dxz_s_s_s_2 = PAx*pz_s_s_s_2 + WPx*pz_s_s_s_3
dyy_s_s_s_2 = PAy*py_s_s_s_2 + WPy*py_s_s_s_3 + 1*oo2p*(s_s_s_s_2 - qoppq*s_s_s_s_3)
dyz_s_s_s_2 = PAy*pz_s_s_s_2 + WPy*pz_s_s_s_3
dzz_s_s_s_2 = PAz*pz_s_s_s_2 + WPz*pz_s_s_s_3 + 1*oo2p*(s_s_s_s_2 - qoppq*s_s_s_s_3)

dxx_s_px_s_0 = QCx*dxx_s_s_s_0 + WQx*dxx_s_s_s_1 + 1*oo2pq*(px_s_s_s_1)
dxx_s_py_s_0 = QCy*dxx_s_s_s_0 + WQy*dxx_s_s_s_1
dxx_s_pz_s_0 = QCz*dxx_s_s_s_0 + WQz*dxx_s_s_s_1
dxy_s_px_s_0 = QCx*dxy_s_s_s_0 + WQx*dxy_s_s_s_1
dxy_s_py_s_0 = QCy*dxy_s_s_s_0 + WQy*dxy_s_s_s_1
dxy_s_pz_s_0 = QCz*dxy_s_s_s_0 + WQz*dxy_s_s_s_1
dxz_s_px_s_0 = QCx*dxz_s_s_s_0 + WQx*dxz_s_s_s_1
dxz_s_py_s_0 = QCy*dxz_s_s_s_0 + WQy*dxz_s_s_s_1
dxz_s_pz_s_0 = QCz*dxz_s_s_s_0 + WQz*dxz_s_s_s_1
dyy_s_px_s_0 = QCx*dyy_s_s_s_0 + WQx*dyy_s_s_s_1
dyy_s_py_s_0 = QCy*dyy_s_s_s_0 + WQy*dyy_s_s_s_1 + 1*oo2pq*(py_s_s_s_1)
dyy_s_pz_s_0 = QCz*dyy_s_s_s_0 + WQz*dyy_s_s_s_1
dyz_s_px_s_0 = QCx*dyz_s_s_s_0 + WQx*dyz_s_s_s_1
dyz_s_py_s_0 = QCy*dyz_s_s_s_0 + WQy*dyz_s_s_s_1
dyz_s_pz_s_0 = QCz*dyz_s_s_s_0 + WQz*dyz_s_s_s_1
dzz_s_px_s_0 = QCx*dzz_s_s_s_0 + WQx*dzz_s_s_s_1
dzz_s_py_s_0 = QCy*dzz_s_s_s_0 + WQy*dzz_s_s_s_1
dzz_s_pz_s_0 = QCz*dzz_s_s_s_0 + WQz*dzz_s_s_s_1 + 1*oo2pq*(pz_s_s_s_1)

dxx_s_px_s_1 = QCx*dxx_s_s_s_1 + WQx*dxx_s_s_s_2 + 1*oo2pq*(px_s_s_s_2)
dxx_s_py_s_1 = QCy*dxx_s_s_s_1 + WQy*dxx_s_s_s_2
dxx_s_pz_s_1 = QCz*dxx_s_s_s_1 + WQz*dxx_s_s_s_2
dxy_s_px_s_1 = QCx*dxy_s_s_s_1 + WQx*dxy_s_s_s_2
dxy_s_py_s_1 = QCy*dxy_s_s_s_1 + WQy*dxy_s_s_s_2
dxy_s_pz_s_1 = QCz*dxy_s_s_s_1 + WQz*dxy_s_s_s_2
dxz_s_px_s_1 = QCx*dxz_s_s_s_1 + WQx*dxz_s_s_s_2
dxz_s_py_s_1 = QCy*dxz_s_s_s_1 + WQy*dxz_s_s_s_2
dxz_s_pz_s_1 = QCz*dxz_s_s_s_1 + WQz*dxz_s_s_s_2
dyy_s_px_s_1 = QCx*dyy_s_s_s_1 + WQx*dyy_s_s_s_2
dyy_s_py_s_1 = QCy*dyy_s_s_s_1 + WQy*dyy_s_s_s_2 + 1*oo2pq*(py_s_s_s_2)
dyy_s_pz_s_1 = QCz*dyy_s_s_s_1 + WQz*dyy_s_s_s_2
dyz_s_px_s_1 = QCx*dyz_s_s_s_1 + WQx*dyz_s_s_s_2
dyz_s_py_s_1 = QCy*dyz_s_s_s_1 + WQy*dyz_s_s_s_2
dyz_s_pz_s_1 = QCz*dyz_s_s_s_1 + WQz*dyz_s_s_s_2
dzz_s_px_s_1 = QCx*dzz_s_s_s_1 + WQx*dzz_s_s_s_2
dzz_s_py_s_1 = QCy*dzz_s_s_s_1 + WQy*dzz_s_s_s_2
dzz_s_pz_s_1 = QCz*dzz_s_s_s_1 + WQz*dzz_s_s_s_2 + 1*oo2pq*(pz_s_s_s_2)

dxx_s_dxx_s_0 = QCx*dxx_s_px_s_0 + WQx*dxx_s_px_s_1 + 1*oo2q*(dxx_s_s_s_0 - qoppq*dxx_s_s_s_1) + 1*oo2pq*(px_s_px_s_1)
dxx_s_dxy_s_0 = QCx*dxx_s_py_s_0 + WQx*dxx_s_py_s_1 + 1*oo2pq*(px_s_py_s_1)
dxx_s_dxz_s_0 = QCx*dxx_s_pz_s_0 + WQx*dxx_s_pz_s_1 + 1*oo2pq*(px_s_pz_s_1)
dxx_s_dyy_s_0 = QCy*dxx_s_py_s_0 + WQy*dxx_s_py_s_1 + 1*oo2q*(dxx_s_s_s_0 - qoppq*dxx_s_s_s_1)
dxx_s_dyz_s_0 = QCy*dxx_s_pz_s_0 + WQy*dxx_s_pz_s_1
dxx_s_dzz_s_0 = QCz*dxx_s_pz_s_0 + WQz*dxx_s_pz_s_1 + 1*oo2q*(dxx_s_s_s_0 - qoppq*dxx_s_s_s_1)
dxy_s_dxx_s_0 = QCx*dxy_s_px_s_0 + WQx*dxy_s_px_s_1 + 1*oo2q*(dxy_s_s_s_0 - qoppq*dxy_s_s_s_1)
dxy_s_dxy_s_0 = QCx*dxy_s_py_s_0 + WQx*dxy_s_py_s_1
dxy_s_dxz_s_0 = QCx*dxy_s_pz_s_0 + WQx*dxy_s_pz_s_1
dxy_s_dyy_s_0 = QCy*dxy_s_py_s_0 + WQy*dxy_s_py_s_1 + 1*oo2q*(dxy_s_s_s_0 - qoppq*dxy_s_s_s_1)
dxy_s_dyz_s_0 = QCy*dxy_s_pz_s_0 + WQy*dxy_s_pz_s_1
dxy_s_dzz_s_0 = QCz*dxy_s_pz_s_0 + WQz*dxy_s_pz_s_1 + 1*oo2q*(dxy_s_s_s_0 - qoppq*dxy_s_s_s_1)
dxz_s_dxx_s_0 = QCx*dxz_s_px_s_0 + WQx*dxz_s_px_s_1 + 1*oo2q*(dxz_s_s_s_0 - qoppq*dxz_s_s_s_1)
dxz_s_dxy_s_0 = QCx*dxz_s_py_s_0 + WQx*dxz_s_py_s_1
dxz_s_dxz_s_0 = QCx*dxz_s_pz_s_0 + WQx*dxz_s_pz_s_1
dxz_s_dyy_s_0 = QCy*dxz_s_py_s_0 + WQy*dxz_s_py_s_1 + 1*oo2q*(dxz_s_s_s_0 - qoppq*dxz_s_s_s_1)
dxz_s_dyz_s_0 = QCy*dxz_s_pz_s_0 + WQy*dxz_s_pz_s_1
dxz_s_dzz_s_0 = QCz*dxz_s_pz_s_0 + WQz*dxz_s_pz_s_1 + 1*oo2q*(dxz_s_s_s_0 - qoppq*dxz_s_s_s_1)
dyy_s_dxx_s_0 = QCx*dyy_s_px_s_0 + WQx*dyy_s_px_s_1 + 1*oo2q*(dyy_s_s_s_0 - qoppq*dyy_s_s_s_1)
dyy_s_dxy_s_0 = QCx*dyy_s_py_s_0 + WQx*dyy_s_py_s_1
dyy_s_dxz_s_0 = QCx*dyy_s_pz_s_0 + WQx*dyy_s_pz_s_1
dyy_s_dyy_s_0 = QCy*dyy_s_py_s_0 + WQy*dyy_s_py_s_1 + 1*oo2q*(dyy_s_s_s_0 - qoppq*dyy_s_s_s_1) + 1*oo2pq*(py_s_py_s_1)
dyy_s_dyz_s_0 = QCy*dyy_s_pz_s_0 + WQy*dyy_s_pz_s_1 + 1*oo2pq*(py_s_pz_s_1)
dyy_s_dzz_s_0 = QCz*dyy_s_pz_s_0 + WQz*dyy_s_pz_s_1 + 1*oo2q*(dyy_s_s_s_0 - qoppq*dyy_s_s_s_1)
dyz_s_dxx_s_0 = QCx*dyz_s_px_s_0 + WQx*dyz_s_px_s_1 + 1*oo2q*(dyz_s_s_s_0 - qoppq*dyz_s_s_s_1)
dyz_s_dxy_s_0 = QCx*dyz_s_py_s_0 + WQx*dyz_s_py_s_1
dyz_s_dxz_s_0 = QCx*dyz_s_pz_s_0 + WQx*dyz_s_pz_s_1
dyz_s_dyy_s_0 = QCy*dyz_s_py_s_0 + WQy*dyz_s_py_s_1 + 1*oo2q*(dyz_s_s_s_0 - qoppq*dyz_s_s_s_1)
dyz_s_dyz_s_0 = QCy*dyz_s_pz_s_0 + WQy*dyz_s_pz_s_1
dyz_s_dzz_s_0 = QCz*dyz_s_pz_s_0 + WQz*dyz_s_pz_s_1 + 1*oo2q*(dyz_s_s_s_0 - qoppq*dyz_s_s_s_1)
dzz_s_dxx_s_0 = QCx*dzz_s_px_s_0 + WQx*dzz_s_px_s_1 + 1*oo2q*(dzz_s_s_s_0 - qoppq*dzz_s_s_s_1)
dzz_s_dxy_s_0 = QCx*dzz_s_py_s_0 + WQx*dzz_s_py_s_1
dzz_s_dxz_s_0 = QCx*dzz_s_pz_s_0 + WQx*dzz_s_pz_s_1
dzz_s_dyy_s_0 = QCy*dzz_s_py_s_0 + WQy*dzz_s_py_s_1 + 1*oo2q*(dzz_s_s_s_0 - qoppq*dzz_s_s_s_1)
dzz_s_dyz_s_0 = QCy*dzz_s_pz_s_0 + WQy*dzz_s_pz_s_1
dzz_s_dzz_s_0 = QCz*dzz_s_pz_s_0 + WQz*dzz_s_pz_s_1 + 1*oo2q*(dzz_s_s_s_0 - qoppq*dzz_s_s_s_1) + 1*oo2pq*(pz_s_pz_s_1)



end
