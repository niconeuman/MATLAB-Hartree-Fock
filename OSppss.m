function gppss = OSppss(RPAValues,RPBValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues)

%PSSS integrals of orders 0 to 3
%Order 0
px_s_s_s_0 = RPAValues(:,1).*gSSSSNValues(:,1)+RWPValues(:,1).*gSSSSNValues(:,2);
py_s_s_s_0 = RPAValues(:,2).*gSSSSNValues(:,1)+RWPValues(:,2).*gSSSSNValues(:,2);
pz_s_s_s_0 = RPAValues(:,3).*gSSSSNValues(:,1)+RWPValues(:,3).*gSSSSNValues(:,2);

%Order 1
px_s_s_s_1 = RPAValues(:,1).*gSSSSNValues(:,2)+RWPValues(:,1).*gSSSSNValues(:,3);
py_s_s_s_1 = RPAValues(:,2).*gSSSSNValues(:,2)+RWPValues(:,2).*gSSSSNValues(:,3);
pz_s_s_s_1 = RPAValues(:,3).*gSSSSNValues(:,2)+RWPValues(:,3).*gSSSSNValues(:,3);

%PPSS integrals of order 0 to 2
%Order 0
px_px_s_s_0 = RPBValues(:,1).*px_s_s_s_0+RWPValues(:,1).*px_s_s_s_1+0.5./pValues.*(gSSSSNValues(:,1)-qValues./ppqValues.*gSSSSNValues(:,2));
px_py_s_s_0 = RPBValues(:,2).*px_s_s_s_0+RWPValues(:,2).*px_s_s_s_1;
px_pz_s_s_0 = RPBValues(:,3).*px_s_s_s_0+RWPValues(:,3).*px_s_s_s_1;

py_px_s_s_0 = RPBValues(:,1).*py_s_s_s_0+RWPValues(:,1).*py_s_s_s_1;
py_py_s_s_0 = RPBValues(:,2).*py_s_s_s_0+RWPValues(:,2).*py_s_s_s_1+0.5./pValues.*(gSSSSNValues(:,1)-qValues./ppqValues.*gSSSSNValues(:,2));
py_pz_s_s_0 = RPBValues(:,3).*py_s_s_s_0+RWPValues(:,3).*py_s_s_s_1;

pz_px_s_s_0 = RPBValues(:,1).*pz_s_s_s_0+RWPValues(:,1).*pz_s_s_s_1;
pz_py_s_s_0 = RPBValues(:,2).*pz_s_s_s_0+RWPValues(:,2).*pz_s_s_s_1;
pz_pz_s_s_0 = RPBValues(:,3).*pz_s_s_s_0+RWPValues(:,3).*pz_s_s_s_1+0.5./pValues.*(gSSSSNValues(:,1)-qValues./ppqValues.*gSSSSNValues(:,2));

gppss = zeros(3,3);
gppss(1,1) = sum(px_px_s_s_0);
gppss(1,2) = sum(px_py_s_s_0);
gppss(1,3) = sum(px_pz_s_s_0);

gppss(2,1) = sum(py_px_s_s_0);
gppss(2,2) = sum(py_py_s_s_0);
gppss(2,3) = sum(py_pz_s_s_0);

gppss(3,1) = sum(pz_px_s_s_0);
gppss(3,2) = sum(pz_py_s_s_0);
gppss(3,3) = sum(pz_pz_s_s_0);

end
