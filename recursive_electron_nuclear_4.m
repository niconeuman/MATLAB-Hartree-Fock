function [VL1L2_N,VL1L2_Np1,bufferInts] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1,L2,Lmax,order,Boys_Table,bufferInts)

%January 6th 2017
%This function calculates the electron nuclear attraction matrices recursively
%The most difficult part is to generate many integrals in a batch.
%For example if I have <p|s> shell doublets, I will want to generate a 3 x
%1 column vector, if I have <p|p> I need a 3 x 3 matrix, for <d|p> a 6 x 3
%(if I use 6 d functions) matrix, etc. How to automatically determine the
%shape and the way to calculate all matrix elements without calling the
%function more times than necessary is the difficult part. But this will be
%much more simple than for the 2 electron integrals, so it is a good time
%to practice.

%Data I need
%L1 = total angular momentum of the first basis function
%L2 = total angular momentum of the second basis function

%According to Helgaker's book, p. 378
Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
Dim2 = (L2+1)*(L2+2)/2;

VL1L2_N = zeros(Dim1,Dim2);
VL1L2_Np1 = zeros(Dim1,Dim2);


if (L1 == 0 && L2 == 0)
    %In each block I need to make the function check if the bufferInts{ind}
    %contains a double. If it does, that will be the output of the
    %function. If it doesn't, it will calculate it.
    %
    %Kab = gprod(RAB(1),RAB(2),RAB(3),a,RPB(1),RPB(2),RPB(3),b);
    q = a*b/p;
    Kab = exp(-q*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
    
    %31 jan 2017. Consider taking out the -1*Z*(2*pi) out to the outside
    %loop in Build_Nuclear_Attraction.
    
%     for ord = 0:Lmax
%     Vss_m = -1*Z*(2*pi/p)*Interpolated_Boys(ord,p*RPC2,Boys_Table)*Kab;
%     %Generally (ind=L1*L2max*(Lmax)+L2*(Lmax)+(m+1))
%     ind = order+1; %As L1 and L2 are zero, I only need the final part of the compound index
%     bufferInts{ind} = Vss_m;
%     end
%     VAB_L1_L2_order = bufferInts{ind};
%     VAB_L1_L2_order_p1 = bufferInts{ind+1};

VL1L2_N = -1*Z*(2*pi/p)*Interpolated_Boys(order,p*RPC2,Boys_Table)*Kab;
VL1L2_Np1 = -1*Z*(2*pi/p)*Interpolated_Boys(order+1,p*RPC2,Boys_Table)*Kab;

elseif (L1 > 0 && L2 == 0)
    nz = [1 4 6 7 9 10 11 13 14 15 16 18 19 20 21 22 24 25 26 27 28 29 31 32 33 34 35 36 37 39 40 41 42 43 44 45 46 48 49 50 51 52 53 54 55 56];
    %This line of integrals will always be a COLUMN vector
    [VL1m1L2_N,VL1m1L2_Np1] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,Lmax,order,Boys_Table,bufferInts);
    [~,VL1m1L2_Np2] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,Lmax,order+1,Boys_Table,bufferInts);
    switch L1-2
        case -1
            
            
            VL1L2_N = [RPA(1)*VL1m1L2_N;RPA(2)*VL1m1L2_N(end-L1+1:end);RPA(3)*VL1m1L2_N(end)]...
                        -[RPC(1)*VL1m1L2_Np1;RPC(2)*VL1m1L2_Np1(end-L1+1:end);RPC(3)*VL1m1L2_Np1(end)];
            VL1L2_Np1 = [RPA(1)*VL1m1L2_Np1;RPA(2)*VL1m1L2_Np1(end-L1+1:end);RPA(3)*VL1m1L2_Np1(end)]...
                        -[RPC(1)*VL1m1L2_Np2;RPC(2)*VL1m1L2_Np2(end-L1+1:end);RPC(3)*VL1m1L2_Np2(end)];
        case 0
            [VL1m2L2_N,VL1m2L2_Np1] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,Lmax,order,Boys_Table,bufferInts);
        %The outputs are [ss](m), and the final integral VAB_L1_L2_order =
        %[ds]
            VL1m2L2_N = VL1m2L2_N*ones(3,1);
            VL1m2L2_Np1 = VL1m2L2_Np1*ones(3,1);
            nonzero = zeros(Dim1,Dim2);
            nonzero(nz(1:length(VL1m2L2_N))) = 1;
            
            Diff_VL1m2L2 = ([VL1m2L2_N-VL1m2L2_Np1;VL1m2L2_N(end-L1+1:end)-VL1m2L2_Np1(end-L1+1:end);VL1m2L2_N(end)-VL1m2L2_Np1(end)]);
            VL1L2_N = [RPA(1)*VL1m1L2_N;RPA(2)*VL1m1L2_N(end-L1+1:end);RPA(3)*VL1m1L2_N(end)]...
                        -[RPC(1)*VL1m1L2_Np1;RPC(2)*VL1m1L2_Np1(end-L1+1:end);RPC(3)*VL1m1L2_Np1(end)]...
                        +1/2/p*(L1-1)*(Diff_VL1m2L2.*nonzero);

        otherwise
            [VL1m2L2_N,VL1m2L2_Np1] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,Lmax,order,Boys_Table,bufferInts);
        %The outputs are [ps](m), and the final integral is [fs] or
        %[ds](m) and the final integrals [gs], etc
        %With this I generate an array the same size as VAB_L1_m1_order(_p1)
            VL1m2L2_N = [VL1m2L2_N;VL1m2L2_N(end-(L1-1)+1:end);VL1m2L2_N(end)];
            
            nonzero = zeros(Dim1,Dim2);
            nonzero(nz(1:length(VL1m2L2_N))) = 1;
            
            Diff_VL1m2L2 = ([VL1m2L2_N-VL1m2L2_Np1;VL1m2L2_N(end-L1+1:end)-VL1m2L2_Np1(end-L1+1:end);VL1m2L2_N(end)-VL1m2L2_Np1(end)]);
            
            VL1L2_N = [RPA(1)*VL1m1L2_N;RPA(2)*VL1m1L2_N(end-L1+1:end);RPA(3)*VL1m1L2_N(end)]...
                        -[RPC(1)*VL1m1L2_Np1;RPC(2)*VL1m1L2_Np1(end-L1+1:end);RPC(3)*VL1m1L2_Np1(end)]...
                        +1/2/p*(L1-1)*(Diff_VL1m2L2.*nonzero);
    end
            

                    
elseif (L1 == 0 && L2 > 0)
    nz = [1 4 6 7 9 10 11 13 14 15 16 18 19 20 21 22 24 25 26 27 28 29 31 32 33 34 35 36 37 39 40 41 42 43 44 45 46 48 49 50 51 52 53 54 55 56]; %later I need to extend this
    %This line of integrals will always be a COLUMN vector
    [VL1L2m1_N,VL1L2m1_Np1] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1,L2-1,Lmax,order,Boys_Table,bufferInts);
    [~,VL1L2m1_Np2] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1,L2-1,Lmax,order+1,Boys_Table,bufferInts);
    switch L2-2
        case -1
            VL1L2_N = [RPB(1)*VL1L2m1_N,RPB(2)*VL1L2m1_N(end-L2+1:end),RPB(3)*VL1L2m1_N(end)]...
                            -[RPC(1)*VL1L2m1_Np1,RPC(2)*VL1L2m1_Np1(end-L2+1:end),RPC(3)*VL1L2m1_Np1(end)];
            VL1L2_Np1 = [RPB(1)*VL1L2m1_Np1,RPB(2)*VL1L2m1_Np1(end-L2+1:end),RPB(3)*VL1L2m1_Np1(end)]...
                            -[RPC(1)*VL1L2m1_Np2,RPC(2)*VL1L2m1_Np2(end-L2+1:end),RPC(3)*VL1L2m1_Np2(end)];
        case 0
            [VL1L2m2_N,VL1L2m2_Np1] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1,L2-2,Lmax,order,Boys_Table,bufferInts);
        %The outputs are [ss](m), and the final integral VAB_L1_L2_order =
        %[ds]
            VL1L2m2_N = VL1L2m2_N*ones(1,3);
            VL1L2m2_Np1 = VL1L2m2_Np1*ones(1,3);
            nonzero = zeros(Dim1,Dim2);
            nonzero(nz(1:length(VL1L2m2_N))) = 1;
            
            Diff_VL1L2m2 = ([VL1L2m2_N-VL1L2m2_Np1,VL1L2m2_N(end-L2+1:end)-VL1L2m2_Np1(end-L2+1:end),VL1L2m2_N(end)-VL1L2m2_Np1(end)]);
            
            VL1L2_N = [RPB(1)*VL1L2m1_N,RPB(2)*VL1L2m1_N(end-L2+1:end),RPB(3)*VL1L2m1_N(end)]...
                            -[RPC(1)*VL1L2m1_Np1,RPC(2)*VL1L2m1_Np1(end-L2+1:end),RPC(3)*VL1L2m1_Np1(end)]...
                            +1/2/p*(L2-1)*(Diff_VL1L2m2.*nonzero);    
        otherwise
            [VL1L2m2_N,VL1L2m2_Np1] = recursive_electron_nuclear_4(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1,L2-2,Lmax,order,Boys_Table,bufferInts);
        %The outputs are [ps](m), and the final integral is [fs] or
        %[ds](m) and the final integrals [gs], etc
        %With this I generate an array the same size as VAB_L1_m1_order(_p1)
            VL1L2m2_N = [VL1L2m2_N,VL1L2m2_N(end-(L2-1)+1:end),VL1L2m2_N(end)]; 
        %15 feb 2017. CAREFUL:
        %If I'm calculating [xx xy xz yy yz zz], (ds), the L2-2 part only adds up
        %to components (1) (4) (6)
        %If I'm calculating [xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz], (fs), the L2-2 part only adds up
        %to components (1) (4) (6) (7) (9) (10)
        %If I'm calculating [xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz], (gs), the L2-2 part only adds up
        %to components (1) (4) (6) (7) (9) (10) (11) (13) (14) (15)
        %If I'm calculating [xxxxx xxxxy xxxxz xxxyy xxxyz xxxzz xxyyy xxyyz xxyzz xxzzz xyyyy xyyyz xyyzz xyzzz xzzzz, yyyyy yyyyz yyyzz yyzzz yzzzz zzzzz], (hs), the L2-2 part only adds up
        %to components (1) (4) (6) (7) (9) (10) (11) (13) (14) (15) (16) (18) (19) (20) (21)
        %If I'm calculating [Xxxxxx Xxxxxy Xxxxxz Xxxxyy Xxxxyz Xxxxzz Xxxyyy Xxxyyz Xxxyzz Xxxzzz Xxyyyy Xxyyyz Xxyyzz Xxyzzz Xxzzzz, Xyyyyy Xyyyyz Xyyyzz Xyyzzz Xyzzzz Xzzzzz, Yyyyyy Yyyyyz Yyyyzz Yyyzzz Yyzzzz Yzzzzz, Zzzzzz], (is), the L2-2 part only adds up
        %to components (1) (4) (6) (7) (9) (10) (11) (13) (14) (15) (16)
        %(18) (19) (20) (21) (22) (24) (25) (26) (27) (28)
        %The number of non-zero components can be obtained from an index
        %vector (for example for gs, which is the maximum I need to
        %calculate (dd) integrals) taken up to number Dim(L-1). i.e. hs
        %goes up to 21, but contains 15 element (length(gs)). gs goes up to
        %15, and contains 10 elements (length(fs))
        %
            nonzero = zeros(Dim1,Dim2);
            nonzero(nz(1:length(VL1L2m2_N))) = 1;
        Diff_VL1L2m2 = ([VL1L2m2_N-VL1L2m2_Np1,VL1L2m2_N(end-L2+1:end)-VL1L2m2_Np1(end-L2+1:end),VL1L2m2_N(end)-VL1L2m2_Np1(end)]);

        VL1L2_N = [RPB(1)*VL1L2m1_N,RPB(2)*VL1L2m1_N(end-L2+1:end),RPB(3)*VL1L2m1_N(end)]...
                            -[RPC(1)*VL1L2m1_Np1,RPC(2)*VL1L2m1_Np1(end-L2+1:end),RPC(3)*VL1L2m1_Np1(end)]...
                            +1/2/p*(L2-1)*(Diff_VL1L2m2.*nonzero);    
    end
end

end