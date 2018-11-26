function gSSSSNValues = primitiveFactorsSSSS_2(basis_a,basis_b,basis_c,basis_d,Boys_Table)

Nintsab = basis_a.n*basis_b.n;
Nintscd = basis_c.n*basis_d.n;
Nints = Nintsab*Nintscd;

gSSSSNValues = zeros(Nints,1);
xValues = zeros(Nints,1);
indexValues = zeros(Nints,1);
xIndexValues = zeros(Nints,1);
DxValues = zeros(Nints,1);

KabValues = zeros(Nintsab,1);
KcdValues = zeros(Nintscd,1);

pValues = zeros(Nintsab,1);
qValues = zeros(Nintscd,1);

WeightValuesab = zeros(Nintsab,1);
WeightValuescd = zeros(Nintscd,1);

PxValues = zeros(Nintsab,1);
PyValues = zeros(Nintsab,1);
PzValues = zeros(Nintsab,1);

QxValues = zeros(Nintscd,1);
QyValues = zeros(Nintscd,1);
QzValues = zeros(Nintscd,1);

alphaValues = zeros(Nints,1);
RPQ2Values = zeros(Nints,1);
PrefactorValues = zeros(Nints,1);




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

        p = aa + ab;
        Px = (aa*g1.x0 + ab*g2.x0)/p;
        Py = (aa*g1.y0 + ab*g2.y0)/p;
        Pz = (aa*g1.z0 + ab*g2.z0)/p;

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
        pValues(t) = p;
        PxValues(t) = Px;
        PyValues(t) = Py;
        PzValues(t) = Pz;
        KabValues(t) = Kab;
        WeightValuesab(t) = c1N1*c2*N2;
        t = t + 1;
    end
end

s = 1;
for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
      g3 = basis_c.g(nc);
      ac = g3.alpha;
      c3 = basis_c.c(nc);
      N3 = g3.N;
      c3N3 = c3*N3;

      for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
            g4 = basis_d.g(nd);
            ad = g4.alpha;
            c4 = basis_d.c(nd);
            N4 = g4.N;

            q = ac+ad;
            Qx = (ac*g3.x0 + ad*g4.x0)/q;
            Qy = (ac*g3.y0 + ad*g4.y0)/q;
            Qz = (ac*g3.z0 + ad*g4.z0)/q;

            rhoCD = ac*ad/q;
            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

            qValues(s) = q;
            QxValues(s) = Qx;
            QyValues(s) = Qy;
            QzValues(s) = Qz;
            KcdValues(s) = Kcd;
            WeightValuescd(s) = c3N3*c4*N4;

            s = s + 1;
      end
end

%I have vectors for all combinations a,b and vectors for all combinations c,d
%for the different quantities I need. %I will construct a matrix which then
%I will cast as a vector of dimension Nints,1

RPQx = PxValues'-QxValues; %' %This will give a matrix
RPQy = PyValues'-QyValues; %' %This will give a matrix
RPQz = PzValues'-QzValues; %' %This will give a matrix

RPQ2 = RPQx.*RPQx+RPQy.*RPQy+RPQz.*RPQz;

ppqValues = qValues + pValues'; %' This should be a matrix
%disp(ppqValues);
pqValues = qValues*pValues'; %'
alpha = pqValues./ppqValues;

x = alpha.*RPQ2;

KabcdValues = KcdValues*KabValues'; %'
%disp(KabcdValues);
PrefactorValues = KabcdValues.*2*17.493418327624862./(pqValues.*sqrt(ppqValues));
PrefactorValues = PrefactorValues(:);
%disp(PrefactorValues);
WeightValues = WeightValuescd*WeightValuesab'; %'
WeightValues = WeightValues(:);
%disp(WeightValues);
xValues = x(:);
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
