function gspss = OSspss(RPBValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues)

%I am changing these functions so that all the primitive shell specific operations are done in a loop in shellOS
%Then the OSabcd functions work in a vectorized manner

%This is a step by step transcription of Table 1 in the Obara-Saika paper

%SPSS integrals of orders 0 to 1
%Order 0
s_px_s_s_0 = RPBValues(:,1).*gSSSSNValues(:,1)+RWPValues(:,1).*gSSSSNValues(:,2);
s_py_s_s_0 = RPBValues(:,2).*gSSSSNValues(:,1)+RWPValues(:,2).*gSSSSNValues(:,2);
s_pz_s_s_0 = RPBValues(:,3).*gSSSSNValues(:,1)+RWPValues(:,3).*gSSSSNValues(:,2);

gpsss = zeros(3,1,1,1);

gspss(1,1,1,1) = sum(s_px_s_s_0);
gspss(1,2,1,1) = sum(s_py_s_s_0);
gspss(1,3,1,1) = sum(s_pz_s_s_0);

end
