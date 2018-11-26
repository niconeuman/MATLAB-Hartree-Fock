function [RPAValues,RPBValues,RPCValues,pValues,VssNValues] = primitiveFactorsNuc(basis_a,basis_b,L1,L2,Boys_Table,AL,Z)

nAtom = size(Z,2);

Nints = basis_a.n*basis_b.n*nAtom;
%boysValues = zeros(Nints,(L1+L2+1));
VssNValues = zeros(Nints,(L1+L2+1));
xValues = zeros(Nints,1);
indexValues = zeros(Nints,1);
xIndexValues = zeros(Nints,1);
DxValues = zeros(Nints,1);
KabValues = zeros(Nints,1);
RPAValues = zeros(Nints,3);
RPBValues = zeros(Nints,3);
RPCValues = zeros(Nints,3);
pValues = zeros(Nints,1);
alphaValues = zeros(Nints,1);
RPQ2Values = zeros(Nints,1);
PrefactorValues = zeros(Nints,1);
WeightValues = zeros(Nints,1);
ZValues = zeros(Nints,1);

xstep = 0.1; %Nlast/Npoints
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
        RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
        RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

        rhoAB = aa*ab/p;
        Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));

        for N = 1:nAtom

                    RPC = [Px-AL(N,1),Py-AL(N,2),Pz-AL(N,3)]; %column vector
                    RPC2 = RPC(1)^2+RPC(2)^2+RPC(3)^2;

                    x = p*RPC2;

                    ZValues(t) = Z(N);
                    xValues(t) = x;
                    KabValues(t) = Kab;
                    RPAValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    RPBValues(t,:) = [RPA(1),RPA(2),RPA(3)];
                    RPCValues(t,:) = [RPC(1),RPC(2),RPC(3)];
                    pValues(t) = p;
                    WeightValues(t) = c1*N1*c2*N2;

                    t = t + 1;


            end
        end
    end

%xValues can be larger than what has been tabulated, so I divide the vector into two conditions

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


PrefactorValues = -KabValues.*ZValues.*6.28318530717959./pValues;

boysValues = zeros(length(xValuesLow),(L1+L2+1));
    for order = 0:(L1+L2)
        boysValues(:,order+1) = Boys_Table(indexValues,order+1)-Boys_Table(indexValues,order+2).*DxValues+...
                                Boys_Table(indexValues,order+3).*Dx2_o_2_Values-Boys_Table(indexValues,order+4).*Dx3_o_6_Values+...
                                Boys_Table(indexValues,order+5).*Dx4_o_24_Values-Boys_Table(indexValues,order+6).*Dx5_o_120_Values;
    end

if ~isempty(xValuesHigh)
    sqrt_pi = 1.772453850905516;

    boysValuesHigh_0 = 0.5*sqrt_pi./sqrt(xValuesHigh); %this calculates the order 0 boys function for x > 35

    %The (prod(n+1:2*n)) part generates a vector and multiplies
    %its components. So to use it on a vector, which should generate a matrix,
    %seems to require for loops. So it is better to precompute it and store it
    %on a vector.
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
    boysValuesHigh_1_N = zeros(length(xValuesHigh),L1+L2);

    for order = 1:(L1+L2)
        boysValuesHigh_1_N(:,order) = Prefactor(order)./xValuesHigh.^(order+0.5);
    end

    boysValuesHigh_0_N = [boysValuesHigh_0,boysValuesHigh_1_N];


    boysValues = [boysValues;boysValuesHigh_0_N];
end

    for order = 0:(L1+L2)
        VssNValues(:,order+1) = WeightValues.*PrefactorValues.*boysValues(:,order+1);
    end



%Outputs: RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,ppqValues,gSSSSNValues

end
