function [RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table)

Nints = basis_a.n*basis_b.n*basis_c.n*basis_d.n;
boysValues = zeros(Nints,(L1+L2+L3+L4+1));
gSSSSNValues = zeros(Nints,(L1+L2+L3+L4+1));
xValues = zeros(Nints,1);
indexValues = zeros(Nints,1);
xIndexValues = zeros(Nints,1);
DxValues = zeros(Nints,1);
KabValues = zeros(Nints,1);
KcdValues = zeros(Nints,1);
RPAValues = zeros(Nints,3);
RPBValues = zeros(Nints,3);
RQCValues = zeros(Nints,3);
RQDValues = zeros(Nints,3);
RWPValues = zeros(Nints,3);
RWQValues = zeros(Nints,3);
pValues = zeros(Nints,1);
qValues = zeros(Nints,1);
alphaValues = zeros(Nints,1);
RPQ2Values = zeros(Nints,1);
PrefactorValues = zeros(Nints,1);
WeightValues = zeros(Nints,1);
RAB = [basis_a.g(1).x0-basis_b.g(1).x0;basis_a.g(1).y0-basis_b.g(1).y0;basis_a.g(1).z0-basis_b.g(1).z0];
RCD = [basis_c.g(1).x0-basis_d.g(1).x0;basis_c.g(1).y0-basis_d.g(1).y0;basis_c.g(1).z0-basis_d.g(1).z0];

xstep = 0.1; %Nlast/Npoints
t = 1;

for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
    g1 = basis_a.g(na);
    aa = g1.alpha;
    c1 = basis_a.c(na);
    N1 = g1.N;
    c1N1 = c1*N1;
    for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
        g2 = basis_b.g(nb);
        ab = g2.alpha;
        c2 = basis_b.c(nb);
        N2 = g2.N;
        c1N1c2N2 = c1N1*c2*N2;

        p = aa + ab;
        Px = (aa*g1.x0 + ab*g2.x0)/p;
        Py = (aa*g1.y0 + ab*g2.y0)/p;
        Pz = (aa*g1.z0 + ab*g2.z0)/p;

        RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
        RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
        RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

            for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                g3 = basis_c.g(nc);
                ac = g3.alpha;
                c3 = basis_c.c(nc);
                N3 = g3.N;
                c1N1c2N2c3N3 = c1N1c2N2*c3*N3;

                for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
                    g4 = basis_d.g(nd);
                    ad = g4.alpha;
                    c4 = basis_d.c(nd);
                    N4 = g4.N;

                    q = ac+ad;
                    Qx = (ac*g3.x0 + ad*g4.x0)/q;
                    Qy = (ac*g3.y0 + ad*g4.y0)/q;
                    Qz = (ac*g3.z0 + ad*g4.z0)/q;

                    RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                    RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                    RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];

                    Wx = (p*Px+q*Qx)/(p+q);
                    Wy = (p*Py+q*Qy)/(p+q);
                    Wz = (p*Pz+q*Qz)/(p+q);

                    RWP = [Wx-Px;Wy-Py;Wz-Pz];
                    RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];

                    rhoCD = ac*ad/q;
                    Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

                    RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                    RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                    alpha = q*p/(q+p);

                    x = alpha*RPQ2;

                    xValues(t) = x;
                    KabValues(t) = Kab;
                    KcdValues(t) = Kcd;
                    RPAValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    RPBValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    RQCValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    RQDValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    RWPValues(t,:) = [RWP(1),RWP(2),RWP(3)];
                    RWQValues(t,:) = [RWQ(1),RWQ(2),RWQ(3)];
                    pValues(t) = p;
                    qValues(t) = q;
                    WeightValues(t) = c1N1c2N2c3N3*c4*N4;

                    t = t + 1;

                end
            end
        end
    end

below35 = (xValues < 35);
xValuesLow = xValues(below35);
xValuesHigh = xValues(~below35);

indexValues = floor(xValuesLow/xstep)+1;
xIndexValues = (indexValues-1)*xstep;
DxValues = (xValuesLow-xIndexValues); %Difference which enters the Taylor expansion

Dx2Values = DxValues.*DxValues;
Dx3Values = Dx2Values.*DxValues;
Dx4Values = Dx2Values.*Dx2Values;
Dx5Values = Dx3Values.*Dx2Values;

Dx2_o_2_Values = .5*Dx2Values;
Dx3_o_6_Values = .166666666666667*Dx3Values;
Dx4_o_24_Values = 4.166666666666666e-2*Dx4Values;
Dx5_o_120_Values = 8.333333333333334e-3*Dx5Values;

ppqValues = pValues + qValues;
PrefactorValues = KabValues.*KcdValues.*2*17.493418327624862./(pValues.*qValues.*sqrt(ppqValues));

boysValuesLow = zeros(length(xValuesLow),L1+L2+L3+L4+1);
for order = 0:(L1+L2+L3+L4)
    boysValuesLow(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                               Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                               Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
end

if ~isempty(xValuesHigh)
    sqrt_pi = 1.772453850905516;
    boysValuesHigh_0 = 0.5*sqrt_pi./sqrt(xValuesHigh); %this calculates the order 0 boys function for x > 35

    %I precompute the quantity (prod(n+1:2*n))/2^(2*n+1)*sqrt_pi
    Prefactor = [443.113462726379e-003
                  664.670194089569e-003
                  1.66167548522392e+000
                  5.81586419828372e+000
                  26.1713888922768e+000
                  143.942638907522e+000
                  935.627152898894e+000
                  7.01720364674171e+003
                  59.6462309973045e+003
                  566.639194474393e+003
                  5.94971154198112e+006
                  68.4216827327829e+006
                  855.271034159787e+006
                  11.5461589611571e+009
                  167.419304936778e+009
                  2.59499922652006e+012
                  42.8174872375810e+012
                  749.306026657668e+012
                  13.8621614931669e+015
                  270.312149116754e+015];
    boysValuesHigh_1_N = zeros(length(xValuesHigh),L1+L2+L3+L4);

    for order = 1:(L1+L2+L3+L4)
        boysValuesHigh_1_N(:,order) = Prefactor(order)./xValuesHigh.^(order+0.5);
    end

    boysValuesHigh_0_N = [boysValuesHigh_0,boysValuesHigh_1_N];

    boysValues(below35==1,:) = boysValuesLow;
    boysValues(below35==0,:) = boysValuesHigh_0_N;
    %disp(boysValues);

    %boysValues = [boysValues;boysValuesHigh_0_N];
else
    boysValues = boysValuesLow;
end

for order = 0:(L1+L2+L3+L4)
    gSSSSNValues(:,order+1) = WeightValues.*PrefactorValues.*boysValues(:,order+1);
end
%Outputs: RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,ppqValues,gSSSSNValues

end
