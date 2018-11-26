function gpsps = OSpsps(RPAValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues)

%I am changing these functions so that all the primitive shell specific operations are done in a loop in shellOS
%Then the OSabcd functions work in a vectorized manner

%This is a step by step transcription of Table 1 in the Obara-Saika paper

%PSSS integrals of orders 0 to 1
%Order 0
px_s_s_s_0 = RPAValues(:,1).*gSSSSNValues(:,1)+RWPValues(:,1).*gSSSSNValues(:,2);
py_s_s_s_0 = RPAValues(:,2).*gSSSSNValues(:,1)+RWPValues(:,2).*gSSSSNValues(:,2);
pz_s_s_s_0 = RPAValues(:,3).*gSSSSNValues(:,1)+RWPValues(:,3).*gSSSSNValues(:,2);

%Order 1
px_s_s_s_1 = RPAValues(:,1).*gSSSSNValues(:,2)+RWPValues(:,1).*gSSSSNValues(:,3);
py_s_s_s_1 = RPAValues(:,2).*gSSSSNValues(:,2)+RWPValues(:,2).*gSSSSNValues(:,3);
pz_s_s_s_1 = RPAValues(:,3).*gSSSSNValues(:,2)+RWPValues(:,3).*gSSSSNValues(:,3);

%PSPS integrals of order 0
%Order 0
px_s_px_s_0 = RQCValues(:,1).*px_s_s_s_0+RWQValues(:,1).*px_s_s_s_1+0.5./ppqValues.*gSSSSNValues(:,2);
px_s_py_s_0 = RQCValues(:,2).*px_s_s_s_0+RWQValues(:,2).*px_s_s_s_1;
px_s_pz_s_0 = RQCValues(:,3).*px_s_s_s_0+RWQValues(:,3).*px_s_s_s_1;

py_s_px_s_0 = RQCValues(:,1).*py_s_s_s_0+RWQValues(:,1).*py_s_s_s_1;
py_s_py_s_0 = RQCValues(:,2).*py_s_s_s_0+RWQValues(:,2).*py_s_s_s_1+0.5./ppqValues.*gSSSSNValues(:,2);
py_s_pz_s_0 = RQCValues(:,3).*py_s_s_s_0+RWQValues(:,3).*py_s_s_s_1;

pz_s_px_s_0 = RQCValues(:,1).*pz_s_s_s_0+RWQValues(:,1).*pz_s_s_s_1;
pz_s_py_s_0 = RQCValues(:,2).*pz_s_s_s_0+RWQValues(:,2).*pz_s_s_s_1;
pz_s_pz_s_0 = RQCValues(:,3).*pz_s_s_s_0+RWQValues(:,3).*pz_s_s_s_1+0.5./ppqValues.*gSSSSNValues(:,2);

gpsps = zeros(3,1,3,1);

gpsps(1,1,1,1) = sum(px_s_px_s_0);
gpsps(1,1,2,1) = sum(px_s_py_s_0);
gpsps(1,1,3,1) = sum(px_s_pz_s_0);

gpsps(2,1,1,1) = sum(py_s_px_s_0);
gpsps(2,1,2,1) = sum(py_s_py_s_0);
gpsps(2,1,3,1) = sum(py_s_pz_s_0);

gpsps(3,1,1,1) = sum(pz_s_px_s_0);
gpsps(3,1,2,1) = sum(pz_s_py_s_0);
gpsps(3,1,3,1) = sum(pz_s_pz_s_0);

end
