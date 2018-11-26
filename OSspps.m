function gspps = OSspps(RPBValues,RQCValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues)

%SPSS integrals of orders 0 to 3
%Order 0
s_px_s_s_0 = RPBValues(:,1).*gSSSSNValues(:,1)+RWPValues(:,1).*gSSSSNValues(:,2);
s_py_s_s_0 = RPBValues(:,2).*gSSSSNValues(:,1)+RWPValues(:,2).*gSSSSNValues(:,2);
s_pz_s_s_0 = RPBValues(:,3).*gSSSSNValues(:,1)+RWPValues(:,3).*gSSSSNValues(:,2);

%Order 1
s_px_s_s_1 = RPBValues(:,1).*gSSSSNValues(:,2)+RWPValues(:,1).*gSSSSNValues(:,3);
s_py_s_s_1 = RPBValues(:,2).*gSSSSNValues(:,2)+RWPValues(:,2).*gSSSSNValues(:,3);
s_pz_s_s_1 = RPBValues(:,3).*gSSSSNValues(:,2)+RWPValues(:,3).*gSSSSNValues(:,3);

%SPPS integarls of order 0 to 1
%Order 0
s_px_px_s_0 = RQCValues(:,1).*s_px_s_s_0+RWQValues(:,1).*s_px_s_s_1+0.5./ppqValues.*gSSSSNValues(:,2);
s_px_py_s_0 = RQCValues(:,2).*s_px_s_s_0+RWQValues(:,2).*s_px_s_s_1;
s_px_pz_s_0 = RQCValues(:,3).*s_px_s_s_0+RWQValues(:,3).*s_px_s_s_1;

s_py_px_s_0 = RQCValues(:,1).*s_py_s_s_0+RWQValues(:,1).*s_py_s_s_1;
s_py_py_s_0 = RQCValues(:,2).*s_py_s_s_0+RWQValues(:,2).*s_py_s_s_1+0.5./ppqValues.*gSSSSNValues(:,2);
s_py_pz_s_0 = RQCValues(:,3).*s_py_s_s_0+RWQValues(:,3).*s_py_s_s_1;

s_pz_px_s_0 = RQCValues(:,1).*s_pz_s_s_0+RWQValues(:,1).*s_pz_s_s_1;
s_pz_py_s_0 = RQCValues(:,2).*s_pz_s_s_0+RWQValues(:,2).*s_pz_s_s_1;
s_pz_pz_s_0 = RQCValues(:,3).*s_pz_s_s_0+RWQValues(:,3).*s_pz_s_s_1+0.5./ppqValues.*gSSSSNValues(:,2);

gspps = zeros(1,3,3);

gspps(1,1,1) = sum(s_px_px_s_0);
gspps(1,1,2) = sum(s_px_py_s_0);
gspps(1,1,3) = sum(s_px_pz_s_0);

gspps(1,2,1) = sum(s_py_px_s_0);
gspps(1,2,2) = sum(s_py_py_s_0);
gspps(1,2,3) = sum(s_py_pz_s_0);

gspps(1,3,1) = sum(s_pz_px_s_0);
gspps(1,3,2) = sum(s_pz_py_s_0);
gspps(1,3,3) = sum(s_pz_pz_s_0);


end
