function [KabcdValues,PrefactorValues,WeightValues,RPQ2Values,gSSSSNValues] = primitiveFactorsSSSS(basis_a,basis_b,basis_c,basis_d,Boys_Table)

Nints = basis_a.n*basis_b.n*basis_c.n*basis_d.n;
%boysValues = zeros(Nints,1);
gSSSSNValues = zeros(Nints,1);
xValues = zeros(Nints,1);
indexValues = zeros(Nints,1);
xIndexValues = zeros(Nints,1);
DxValues = zeros(Nints,1);
KabValues = zeros(Nints,1);
KcdValues = zeros(Nints,1);
%RPAValues = zeros(Nints,3);
%RPBValues = zeros(Nints,3);
%RQCValues = zeros(Nints,3);
%RQDValues = zeros(Nints,3);
%RWPValues = zeros(Nints,3);
%RWQValues = zeros(Nints,3);
RPQ2Values = zeros(Nints,3);
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
        %RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

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

                    %RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                    %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                    %RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];

                    %Wx = (p*Px+q*Qx)/(p+q);
                    %Wy = (p*Py+q*Qy)/(p+q);
                    %Wz = (p*Pz+q*Qz)/(p+q);

                    %RWP = [Wx-Px;Wy-Py;Wz-Pz];
                    %RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];

                    rhoCD = ac*ad/q;
                    Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

                    RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                    RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                    RPQ2Values(t) = RPQ2;
                    alpha = q*p/(q+p);

                    x = alpha*RPQ2;

                    xValues(t) = x;
                    KabValues(t) = Kab;
                    KcdValues(t) = Kcd;
                    %RPAValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    %RPBValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    %RQCValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    %RQDValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    %RWPValues(t,:) = [RWP(1),RWP(2),RWP(3)];
                    %RWQValues(t,:) = [RWQ(1),RWQ(2),RWQ(3)];
                    pValues(t) = p;
                    qValues(t) = q;
                    WeightValues(t) = c1N1c2N2c3N3*c4*N4;

                    t = t + 1;

                end
            end
        end
    end
%disp(xValues);

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
%disp(ppqValues);
PrefactorValues = KabValues.*KcdValues.*2*17.493418327624862./(pValues.*qValues.*sqrt(ppqValues));
KabcdValues = KabValues.*KcdValues;
%disp(KabValues.*KcdValues);
%disp(PrefactorValues);
%disp(WeightValues);
boysValues = zeros(length(xValuesLow),1);
order = 0;
    boysValues(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                            Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                            Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;

if ~isempty(xValuesHigh)
    sqrt_pi = 1.772453850905516;

    boysValuesHigh_0 = 0.5*sqrt_pi./sqrt(xValuesHigh); %this calculates the order 0 boys function for x > 35

    boysValues = [boysValues;boysValuesHigh_0];
end

gSSSSNValues = WeightValues.*PrefactorValues.*boysValues;

%Outputs: gSSSSNValues

end
