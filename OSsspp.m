function gsspp = OSsspp(RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues)

%SSPS integrals of orders 0 to 3
%Order 0
s_s_px_s_0 = RQCValues(:,1).*gSSSSNValues(:,1)+RWQValues(:,1).*gSSSSNValues(:,2);
s_s_py_s_0 = RQCValues(:,2).*gSSSSNValues(:,1)+RWQValues(:,2).*gSSSSNValues(:,2);
s_s_pz_s_0 = RQCValues(:,3).*gSSSSNValues(:,1)+RWQValues(:,3).*gSSSSNValues(:,2);

%Order 1
s_s_px_s_1 = RQCValues(:,1).*gSSSSNValues(:,2)+RWQValues(:,1).*gSSSSNValues(:,3);
s_s_py_s_1 = RQCValues(:,2).*gSSSSNValues(:,2)+RWQValues(:,2).*gSSSSNValues(:,3);
s_s_pz_s_1 = RQCValues(:,3).*gSSSSNValues(:,2)+RWQValues(:,3).*gSSSSNValues(:,3);

%PPSS integrals of order 0 to 2
%Order 0
s_s_px_px_0 = RQDValues(:,1).*s_s_px_s_0+RWQValues(:,1).*s_s_px_s_1+0.5./qValues.*(gSSSSNValues(:,1)-pValues./ppqValues.*gSSSSNValues(:,2));
s_s_px_py_0 = RQDValues(:,2).*s_s_px_s_0+RWQValues(:,2).*s_s_px_s_1;
s_s_px_pz_0 = RQDValues(:,3).*s_s_px_s_0+RWQValues(:,3).*s_s_px_s_1;

s_s_py_px_0 = RQDValues(:,1).*s_s_py_s_0+RWQValues(:,1).*s_s_py_s_1;
s_s_py_py_0 = RQDValues(:,2).*s_s_py_s_0+RWQValues(:,2).*s_s_py_s_1+0.5./qValues.*(gSSSSNValues(:,1)-pValues./ppqValues.*gSSSSNValues(:,2));
s_s_py_pz_0 = RQDValues(:,3).*s_s_py_s_0+RWQValues(:,3).*s_s_py_s_1;

s_s_pz_px_0 = RQDValues(:,1).*s_s_pz_s_0+RWQValues(:,1).*s_s_pz_s_1;
s_s_pz_py_0 = RQDValues(:,2).*s_s_pz_s_0+RWQValues(:,2).*s_s_pz_s_1;
s_s_pz_pz_0 = RQDValues(:,3).*s_s_pz_s_0+RWQValues(:,3).*s_s_pz_s_1+0.5./qValues.*(gSSSSNValues(:,1)-pValues./ppqValues.*gSSSSNValues(:,2));

gsspp = zeros(1,1,3,3);

gsspp(1,1,1,1) = sum(s_s_px_px_0);
gsspp(1,1,1,2) = sum(s_s_px_py_0);
gsspp(1,1,1,3) = sum(s_s_px_pz_0);

gsspp(1,1,2,1) = sum(s_s_py_px_0);
gsspp(1,1,2,2) = sum(s_s_py_py_0);
gsspp(1,1,2,3) = sum(s_s_py_pz_0);

gsspp(1,1,3,1) = sum(s_s_pz_px_0);
gsspp(1,1,3,2) = sum(s_s_pz_py_0);
gsspp(1,1,3,3) = sum(s_s_pz_pz_0);

end
