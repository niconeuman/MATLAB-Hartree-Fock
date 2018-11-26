function gabcd = shellOS(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table)

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

xstep = 0.1; %Nlast/Npoints

if ( L1 == 1 && L2 == 0 && L3 == 0 && L4 == 0) %px_s_s_s_0
%tic;
t = 1;
for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
    g1 = basis_a.g(na);
    aa = g1.alpha;
    c1 = basis_a.c(na);
    N1 = g1.N;
    for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
        g2 = basis_b.g(nb);
        ab = g2.alpha;
        c2 = basis_b.c(nb);
        N2 = g2.N;

        p = aa + ab;
        Px = (aa*g1.x0 + ab*g2.x0)/p;
        Py = (aa*g1.y0 + ab*g2.y0)/p;
        Pz = (aa*g1.z0 + ab*g2.z0)/p;

        RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
        %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
        RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

            for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                g3 = basis_c.g(nc);
                ac = g3.alpha;
                c3 = basis_c.c(nc);
                N3 = g3.N;

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
                    RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];

                    Wx = (p*Px+q*Qx)/(p+q);
                    Wy = (p*Py+q*Qy)/(p+q);
                    Wz = (p*Pz+q*Qz)/(p+q);

                    RWP = [Wx-Px;Wy-Py;Wz-Pz];
                    %RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];

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
                    %RQCValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    RWPValues(t,:) = [RWP(1),RWP(2),RWP(3)];
                    %RWQValues(t,:) = [RWQ(1),RWQ(2),RWQ(3)];
                    pValues(t) = p;
                    qValues(t) = q;
                    WeightValues(t) = c1*N1*c2*N2*c3*N3*c4*N4;

                    t = t + 1;

                end
            end
        end
    end

indexValues = floor(xValues/xstep)+1;
xIndexValues = (indexValues-1)*xstep;
DxValues = (xValues-xIndexValues); %Difference which enters the Taylor expansion

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

    for order = 0:(L1+L2+L3+L4)
        boysValues(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                                Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
        gSSSSNValues(:,order+1) = WeightValues.*PrefactorValues.*boysValues(:,order+1);
    end


    gpsss = OSpsss(RPAValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues);
    gabcd = gpsss;
%    disp([basis_a.g(1).x0,basis_b.g(1).x0,basis_c.g(1).y0,basis_d.g(1).x0]);
%    disp([basis_a.g(1).alpha,basis_b.g(1).alpha,basis_c.g(1).alpha,basis_d.g(1).alpha]);
%    disp('WeightValues,PrefactorValues,gSSSSNValues');
%    disp([WeightValues,PrefactorValues,gSSSSNValues]);
%disp('Time for gpsss')
%toc;

elseif ( L1 == 0 && L2 == 1 && L3 == 0 && L4 == 0) %s_px_s_s_0
%tic;
t = 1;
for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
    g1 = basis_a.g(na);
    aa = g1.alpha;
    c1 = basis_a.c(na);
    N1 = g1.N;
    for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
        g2 = basis_b.g(nb);
        ab = g2.alpha;
        c2 = basis_b.c(nb);
        N2 = g2.N;

        p = aa + ab;
        Px = (aa*g1.x0 + ab*g2.x0)/p;
        Py = (aa*g1.y0 + ab*g2.y0)/p;
        Pz = (aa*g1.z0 + ab*g2.z0)/p;

        %RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
        RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
        RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

            for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                g3 = basis_c.g(nc);
                ac = g3.alpha;
                c3 = basis_c.c(nc);
                N3 = g3.N;

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
                    RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];

                    Wx = (p*Px+q*Qx)/(p+q);
                    Wy = (p*Py+q*Qy)/(p+q);
                    Wz = (p*Pz+q*Qz)/(p+q);

                    RWP = [Wx-Px;Wy-Py;Wz-Pz];
                    %RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];

                    rhoCD = ac*ad/q;
                    Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

                    RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                    RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                    alpha = q*p/(q+p);

                    x = alpha*RPQ2;

                    xValues(t) = x;
                    KabValues(t) = Kab;
                    KcdValues(t) = Kcd;
                    %RPAValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    RPBValues(t,:) = [RPB(1),RPB(2),RPB(3)];
                    %RQCValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    RWPValues(t,:) = [RWP(1),RWP(2),RWP(3)];
                    %RWQValues(t,:) = [RWQ(1),RWQ(2),RWQ(3)];
                    pValues(t) = p;
                    qValues(t) = q;
                    WeightValues(t) = c1*N1*c2*N2*c3*N3*c4*N4;

                    t = t + 1;

                end
            end
        end
    end

indexValues = floor(xValues/xstep)+1;
xIndexValues = (indexValues-1)*xstep;
DxValues = (xValues-xIndexValues); %Difference which enters the Taylor expansion

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

    for order = 0:(L1+L2+L3+L4)
        boysValues(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                                Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
        gSSSSNValues(:,order+1) = WeightValues.*PrefactorValues.*boysValues(:,order+1);
    end


    gspss = OSspss(RPBValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues);
    gabcd = gspss;
%    disp([basis_a.g(1).x0,basis_b.g(1).x0,basis_c.g(1).y0,basis_d.g(1).x0]);
%    disp([basis_a.g(1).alpha,basis_b.g(1).alpha,basis_c.g(1).alpha,basis_d.g(1).alpha]);
%    disp('WeightValues,PrefactorValues,gSSSSNValues');
%    disp([WeightValues,PrefactorValues,gSSSSNValues]);
%disp('Time for gspss')
%toc;

elseif ( L1 == 1 && L2 == 0 && L3 == 1 && L4 == 0) %px_s_px_s_0
%tic;
t = 1;
for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
    g1 = basis_a.g(na);
    aa = g1.alpha;
    c1 = basis_a.c(na);
    N1 = g1.N;
    for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
        g2 = basis_b.g(nb);
        ab = g2.alpha;
        c2 = basis_b.c(nb);
        N2 = g2.N;

        p = aa + ab;
        Px = (aa*g1.x0 + ab*g2.x0)/p;
        Py = (aa*g1.y0 + ab*g2.y0)/p;
        Pz = (aa*g1.z0 + ab*g2.z0)/p;

        RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
        %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
        RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

            for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                g3 = basis_c.g(nc);
                ac = g3.alpha;
                c3 = basis_c.c(nc);
                N3 = g3.N;

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
                    %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
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
                    RQCValues(t,:) = [RQC(1),RQC(2),RQC(3)];
                    RWPValues(t,:) = [RWP(1),RWP(2),RWP(3)];
                    RWQValues(t,:) = [RWQ(1),RWQ(2),RWQ(3)];
                    pValues(t) = p;
                    qValues(t) = q;
                    WeightValues(t) = c1*N1*c2*N2*c3*N3*c4*N4;

                    t = t + 1;

                end
            end
        end
    end

indexValues = floor(xValues/xstep)+1;
xIndexValues = (indexValues-1)*xstep;
DxValues = (xValues-xIndexValues); %Difference which enters the Taylor expansion

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

    for order = 0:(L1+L2+L3+L4)
        boysValues(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                                Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
        gSSSSNValues(:,order+1) = WeightValues.*PrefactorValues.*boysValues(:,order+1);
    end


    gpsps = OSpsps(RPAValues,RQCValues,RWPValues,RWQValues,ppqValues,gSSSSNValues);
    gabcd = gpsps;
%    disp([basis_a.g(1).x0,basis_b.g(1).x0,basis_c.g(1).y0,basis_d.g(1).x0]);
%    disp([basis_a.g(1).alpha,basis_b.g(1).alpha,basis_c.g(1).alpha,basis_d.g(1).alpha]);
%    disp('WeightValues,PrefactorValues,gSSSSNValues');
%    disp([WeightValues,PrefactorValues,gSSSSNValues]);
%disp('Time for gpsps')
%toc;
elseif ( L1 == 1 && L2 == 0 && L3 == 1 && L4 == 1)

    gpspp_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table);
    gabcd = permute(gpspp_swap,[3 4 1 2]);

elseif ( L1 == 0 && L2 == 1 && L3 == 1 && L4 == 1)

    gsppp_swap = shellOS(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table);
    gabcd = permute(gsppp_swap,[3 4 1 2]);

elseif( L1 == 1 && L2 == 1 && L3 == 0 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);

gppss = OSppss(RPAValues,RPBValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gppss;

elseif( L1 == 0 && L2 == 0 && L3 == 1 && L4 == 1)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
gsspp = OSsspp(RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gsspp;

elseif( L1 == 0 && L2 == 1 && L3 == 1 && L4 == 0)

[RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues] = primitiveFactors(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table);
gspps = OSspps(RPBValues,RQCValues,RWPValues,RWQValues,pValues,qValues,ppqValues,gSSSSNValues);
gabcd = gspps;

elseif( L1 == 0 && L2 == 1 && L3 == 0 && L4 == 1)

gspsp = zeros(1,3,1,3);

%fun is a function handle. I could use VRR, HRR, etc


    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        gspsp_temp = zeros(1,3,1,3);

            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;

                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;

                %RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
                RPA = -ab/p*RAB;
                RPB = aa/p*RAB;


                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

                gspsp_temp2 = zeros(1,3,1,3);


                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;

                        gspsp_temp3 = zeros(1,3,1,3);

                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.

                            g4 = basis_d.g(nd);
                            ad = g4.alpha;
                            c4 = basis_d.c(nd);
                            N4 = g4.N;

                            q = ac+ad;
                            Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            Qz = (ac*g3.z0 + ad*g4.z0)/q;

%                             A_CD = EcdX*EcdY*EcdZ*basis_c.c(nc)*basis_d.c(nd)*basis_c.g(nc).N*basis_d.g(nd).N;
                            %RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
                            RQC = -ad/q*RCD;
                            RQD = ac/q*RCD;
%                             RCDx = g3.x0-g4.x0;
%                             RCDy = g3.y0-g4.y0;
%                             RCDz = g3.z0-g4.z0;
%                             Wx = (p*Px+q*Qx)/(p+q);
%                             Wy = (p*Py+q*Qy)/(p+q);
%                             Wz = (p*Pz+q*Qz)/(p+q);
                            %Wx - Px = (p*Px+q*Qx)/(p+q)-Px =
                            %(p-(p+q))*Px/(p+q)+q/(p+q)*Qx =
                            %-q/(p+q)*Px+q/(p+q)*Qx = -q/(p+q)(Px-Qx)



                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));



                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RWP = -q/(p+q)*RPQ;
                            RWQ = p/(p+q)*RPQ;
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);

                            %gabcd(a,b,c,d) = gabcd(a,b,c,d)+A_AB*A_CD*Boys(0,alpha*RPQ2)*2*pi^2.5/(p*pp*sqrt(p+pp));
                            %temp3 = temp3 + vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table)*c4*N4;
                            %For checking the different recursion
                            %strategies
                            %temp3 = temp3 + vrr_1old(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,[0;0;0],[0;0;0],[0;0;0],RWP,[0;0;0],p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,L1,order,Boys_Table,nz)*c4*N4;
                            %temp3 = temp3 + vrr_1(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)*c4*N4;

                            gspsp_temp4 = OSspsp(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table);
                            gspsp_temp3 = gspsp_temp3 + gspsp_temp4*c4*N4;

                        end
                        gspsp_temp2 = gspsp_temp2 + gspsp_temp3*c3*N3;

                    end
                    gspsp_temp = gspsp_temp + gspsp_temp2*c2*N2;

            end
            gspsp = gspsp + gspsp_temp*c1*N1;


    end
gabcd = gspsp;

elseif( L1 == 1 && L2 == 0 && L3 == 0 && L4 == 1)

gpssp = zeros(3,1,1,3);

%fun is a function handle. I could use VRR, HRR, etc


    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        gpssp_temp = zeros(3,1,1,3);

            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;

                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;

                %RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
                RPA = -ab/p*RAB;
                RPB = aa/p*RAB;


                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

                gpssp_temp2 = zeros(3,1,1,3);


                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;

                        gpssp_temp3 = zeros(3,1,1,3);

                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.

                            g4 = basis_d.g(nd);
                            ad = g4.alpha;
                            c4 = basis_d.c(nd);
                            N4 = g4.N;

                            q = ac+ad;
                            Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            Qz = (ac*g3.z0 + ad*g4.z0)/q;

%                             A_CD = EcdX*EcdY*EcdZ*basis_c.c(nc)*basis_d.c(nd)*basis_c.g(nc).N*basis_d.g(nd).N;
                            %RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
                            RQC = -ad/q*RCD;
                            RQD = ac/q*RCD;
%                             RCDx = g3.x0-g4.x0;
%                             RCDy = g3.y0-g4.y0;
%                             RCDz = g3.z0-g4.z0;
%                             Wx = (p*Px+q*Qx)/(p+q);
%                             Wy = (p*Py+q*Qy)/(p+q);
%                             Wz = (p*Pz+q*Qz)/(p+q);
                            %Wx - Px = (p*Px+q*Qx)/(p+q)-Px =
                            %(p-(p+q))*Px/(p+q)+q/(p+q)*Qx =
                            %-q/(p+q)*Px+q/(p+q)*Qx = -q/(p+q)(Px-Qx)



                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));



                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RWP = -q/(p+q)*RPQ;
                            RWQ = p/(p+q)*RPQ;
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);

                            %gabcd(a,b,c,d) = gabcd(a,b,c,d)+A_AB*A_CD*Boys(0,alpha*RPQ2)*2*pi^2.5/(p*pp*sqrt(p+pp));
                            %temp3 = temp3 + vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table)*c4*N4;
                            %For checking the different recursion
                            %strategies
                            %temp3 = temp3 + vrr_1old(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,[0;0;0],[0;0;0],[0;0;0],RWP,[0;0;0],p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,L1,order,Boys_Table,nz)*c4*N4;
                            %temp3 = temp3 + vrr_1(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)*c4*N4;

                            gpssp_temp4 = OSpssp(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table);
                            gpssp_temp3 = gpssp_temp3 + gpssp_temp4*c4*N4;

                        end
                        gpssp_temp2 = gpssp_temp2 + gpssp_temp3*c3*N3;

                    end
                    gpssp_temp = gpssp_temp + gpssp_temp2*c2*N2;

            end
            gpssp = gpssp + gpssp_temp*c1*N1;


    end
gabcd = gpssp;

elseif( L1 == 1 && L2 == 1 && L3 == 1 && L4 == 0)

gppps = zeros(3,3,3,1);

    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        gppps_temp = zeros(3,3,3,1);

            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;

                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;

                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
                RPA = -ab/p*RAB;
                RPB = aa/p*RAB;

                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

                gppps_temp2 = zeros(3,3,3,1);

                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;

                        gppps_temp3 = zeros(3,3,3,1);

                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.

                            g4 = basis_d.g(nd);
                            ad = g4.alpha;
                            c4 = basis_d.c(nd);
                            N4 = g4.N;

                            q = ac+ad;
                            Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            Qz = (ac*g3.z0 + ad*g4.z0)/q;

                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
                            RQC = -ad/q*RCD;
                            RQD = ac/q*RCD;

                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RWP = -q/(p+q)*RPQ;
                            RWQ = p/(p+q)*RPQ;
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);

                            gppps_temp4 = OSppps(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table);
                            gppps_temp3 = gppps_temp3 + gppps_temp4*c4*N4;

                        end
                        gppps_temp2 = gppps_temp2 + gppps_temp3*c3*N3;
                    end
                    gppps_temp = gppps_temp + gppps_temp2*c2*N2;
            end
            gppps = gppps + gppps_temp*c1*N1;
    end
gabcd = gppps;

elseif( L1 == 1 && L2 == 1 && L3 == 0 && L4 == 1)

gppsp = zeros(3,3,1,3);

    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        gppsp_temp = zeros(3,3,1,3);

            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;

                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;

                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
                RPA = -ab/p*RAB;
                RPB = aa/p*RAB;

                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

                gppsp_temp2 = zeros(3,3,1,3);

                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;

                        gppsp_temp3 = zeros(3,3,1,3);

                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.

                            g4 = basis_d.g(nd);
                            ad = g4.alpha;
                            c4 = basis_d.c(nd);
                            N4 = g4.N;

                            q = ac+ad;
                            Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            Qz = (ac*g3.z0 + ad*g4.z0)/q;

                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
                            RQC = -ad/q*RCD;
                            RQD = ac/q*RCD;

                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RWP = -q/(p+q)*RPQ;
                            RWQ = p/(p+q)*RPQ;
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);

                            gppsp_temp4 = OSppsp(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table);
                            gppsp_temp3 = gppsp_temp3 + gppsp_temp4*c4*N4;

                        end
                        gppsp_temp2 = gppsp_temp2 + gppsp_temp3*c3*N3;
                    end
                    gppsp_temp = gppsp_temp + gppsp_temp2*c2*N2;
            end
            gppsp = gppsp + gppsp_temp*c1*N1;
    end
gabcd = gppsp;

elseif( L1 == 1 && L2 == 1 && L3 == 1 && L4 == 1)

gpppp = zeros(3,3,3,3);

%fun is a function handle. I could use VRR, HRR, etc


    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        gpppp_temp = zeros(3,3,3,3);

            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;

                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;

                %RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
                RPA = -ab/p*RAB;
                RPB = aa/p*RAB;


                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

                gpppp_temp2 = zeros(3,3,3,3);


                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;

                        gpppp_temp3 = zeros(3,3,3,3);

                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.

                            g4 = basis_d.g(nd);
                            ad = g4.alpha;
                            c4 = basis_d.c(nd);
                            N4 = g4.N;

                            q = ac+ad;
                            Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            Qz = (ac*g3.z0 + ad*g4.z0)/q;

%                             A_CD = EcdX*EcdY*EcdZ*basis_c.c(nc)*basis_d.c(nd)*basis_c.g(nc).N*basis_d.g(nd).N;
                            %RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
                            RQC = -ad/q*RCD;
                            RQD = ac/q*RCD;
%                             RCDx = g3.x0-g4.x0;
%                             RCDy = g3.y0-g4.y0;
%                             RCDz = g3.z0-g4.z0;
%                             Wx = (p*Px+q*Qx)/(p+q);
%                             Wy = (p*Py+q*Qy)/(p+q);
%                             Wz = (p*Pz+q*Qz)/(p+q);
                            %Wx - Px = (p*Px+q*Qx)/(p+q)-Px =
                            %(p-(p+q))*Px/(p+q)+q/(p+q)*Qx =
                            %-q/(p+q)*Px+q/(p+q)*Qx = -q/(p+q)(Px-Qx)



                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));



                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RWP = -q/(p+q)*RPQ;
                            RWQ = p/(p+q)*RPQ;
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);

                            %gabcd(a,b,c,d) = gabcd(a,b,c,d)+A_AB*A_CD*Boys(0,alpha*RPQ2)*2*pi^2.5/(p*pp*sqrt(p+pp));
                            %temp3 = temp3 + vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table)*c4*N4;
                            %For checking the different recursion
                            %strategies
                            %temp3 = temp3 + vrr_1old(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,[0;0;0],[0;0;0],[0;0;0],RWP,[0;0;0],p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,L1,order,Boys_Table,nz)*c4*N4;
                            %temp3 = temp3 + vrr_1(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)*c4*N4;

                            gpppp_temp4 = OSpppp(Kab,Kcd,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,alpha,RPQ2,Boys_Table);
                            gpppp_temp3 = gpppp_temp3 + gpppp_temp4*c4*N4;

                        end
                        gpppp_temp2 = gpppp_temp2 + gpppp_temp3*c3*N3;

                    end
                    gpppp_temp = gpppp_temp + gpppp_temp2*c2*N2;

            end
            gpppp = gpppp + gpppp_temp*c1*N1;


    end

gabcd = gpppp;


end  %if

end
