function gpsss = OSpsss(RPAValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues)

%I am changing these functions so that all the primitive shell specific operations are done in a loop in shellOS
%Then the OSabcd functions work in a vectorized manner

%This is a step by step transcription of Table 1 in the Obara-Saika paper

%PSSS integrals of orders 0 to 1
%Order 0
px_s_s_s_0 = RPAValues(:,1).*gSSSSNValues(:,1)+RWPValues(:,1).*gSSSSNValues(:,2);
py_s_s_s_0 = RPAValues(:,2).*gSSSSNValues(:,1)+RWPValues(:,2).*gSSSSNValues(:,2);
pz_s_s_s_0 = RPAValues(:,3).*gSSSSNValues(:,1)+RWPValues(:,3).*gSSSSNValues(:,2);

gpsss = zeros(3,1,1,1);

gpsss(1,1,1,1) = sum(px_s_s_s_0);
gpsss(2,1,1,1) = sum(py_s_s_s_0);
gpsss(3,1,1,1) = sum(pz_s_s_s_0);

end
