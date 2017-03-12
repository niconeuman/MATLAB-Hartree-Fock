function [VAB_L1_L2_order,VAB_L1_L2_order_p1] = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1,L2,order,Boys_Table)

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

VAB_L1_L2_order = zeros(Dim1,Dim2);
VAB_L1_L2_order_p1 = zeros(Dim1,Dim2);


if (L1 == 0 && L2 == 0)
%     Ex = gprod_1D(RPA(1),a,RPB(1),b);
%     Ey = gprod_1D(RPA(2),a,RPB(2),b);
%     Ez = gprod_1D(RPA(3),a,RPB(3),b);
    Kab = gprod(RPA(1),RPA(2),RPA(3),a,RPB(1),RPB(2),RPB(3),b);
    %Kab = Ex*Ey*Ez;
    VAB_L1_L2_order = -1*Z*(2*pi/p)*Interpolated_Boys(order,p*RPC2,Boys_Table)*Kab;
    VAB_L1_L2_order_p1 = -1*Z*(2*pi/p)*Interpolated_Boys(order+1,p*RPC2,Boys_Table)*Kab;
    
%     VAB_L1_L2_order = -1*Z*(2*pi/p)*Boys(order,p*RPC2)*Kab;
%     VAB_L1_L2_order_p1 = -1*Z*(2*pi/p)*Boys(order+1,p*RPC2)*Kab;
    
elseif (L1 > 0 && L2 == 0)
    %This line of integrals will always be a COLUMN vector
    [VAB_L1_m1_order,VAB_L1_m1_order_p1] = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order,Boys_Table);
    
    if ((L1 - 2) >= 0)
    [VAB_L1_m2_order,VAB_L1_m2_order_p1] = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,order,Boys_Table);
   
    else
    VAB_L1_m2_order = 0;
    VAB_L1_m2_order_p1 = 0;
    end
    %Check this formula carefully because Eq. 9.10.8 in Helgaker's book
    
    % I NEED to transpose the column vectors into rows
    VAB_matrix = RPA'*VAB_L1_m1_order'+1/2/p*(L1-1)*VAB_L1_m2_order'...
        -RPC'*VAB_L1_m1_order_p1-1/2/p*(L1-1)*VAB_L1_m2_order_p1';
    
    if (size(VAB_matrix,2) == 1)
            VAB_L1_L2_order = VAB_matrix;
    elseif (size(VAB_matrix,2) == 3)
        %[xx xy xz;
        % yx yy yz;
        % zx zy zz];
            VAB_L1_L2_order = [VAB_matrix(1,1:3)';VAB_matrix(2,2:3)';VAB_matrix(3,3)];
            %VAB = [xx
            %       xy
            %       xz
            %       yy
            %       yz
            %       zz];
    elseif (size(VAB_matrix,2) == 6) %For now (6/01/2017) I won't use more than f functions, so I don't think I need more for now
        %[Xxx,Xxy,Xxz,Xyy,Xyz,Xzz;
        % Yxx,Yxy,Yxz,Yyy,Yyz,Yzz;
        % Zxx,Zxy,Zxz,Zyy,Zyz,Zzz];
        VAB_L1_L2_order = [VAB_matrix(1,1:6)';VAB_matrix(2,4:6)';VAB_matrix(3,6)];
%     elseif (size(VAB_matrix,2) == 10)
%         %[Xxxx,Xxxy,Xxxz,Xxyy,Xxyz,Xxzz,Xyyy,Xyyz,Xyzz,Xzzz;
%         % Yxxx,Yxxy,Yxxz,Yxyy,Yxyz,Yxzz,Yyyy,Yyyz,Yyzz,Yzzz;
%         % Zxxx,Zxxy,Zxxz,Zxyy,Zxyz,Zxzz,Zyyy,Zyyz,Zyzz,Zzzz];
%         VAB = [VAB_matrix(1,1:10)';VAB_matrix(2,7:10)';VAB_matrix(3,10)]'; %L = 4, (L+1)*(L+2)/2 = 5*6/2 = 15
    end
    
elseif (L1 == 0 && L2 > 0)

        %This line of integrals will always be a ROW vector
    VAB_L2_minus_1_order = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1,L2-1,order,Boys_Table);
    VAB_L2_minus_1_order_plus_1 = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1,L2-1,order+1,Boys_Table);
    if ((L2 - 2) >= 0)
    VAB_L2_minus_2_order = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1,L2-2,order,Boys_Table);
    VAB_L2_minus_2_order_plus_1 = recursive_electron_nuclear_2(a,b,RPA,RPB,RPC,p,RPC2,Z,L1,L2-2,order+1,Boys_Table);
    else
    VAB_L2_minus_2_order = 0;
    VAB_L2_minus_2_order_plus_1 = 0;
    end
    %Check this formula carefully because Eq. 9.10.8 in Helgaker's book
    
    %I DON'T NEED to transpose the integral vectors.
    VAB_matrix = RPB'*VAB_L2_minus_1_order+1/2/p*(L2-1)*VAB_L2_minus_2_order...
        -RPC'*VAB_L2_minus_1_order_plus_1-1/2/p*(L2-1)*VAB_L2_minus_2_order_plus_1;
    
    if (size(VAB_matrix,2) == 1)
            VAB_L1_L2_order = VAB_matrix'; %L = 1
    elseif (size(VAB_matrix,2) == 3)
        %[xx xy xz;
        % yx yy yz;
        % zx zy zz];
            VAB_L1_L2_order = [VAB_matrix(1,1:3)';VAB_matrix(2,2:3)';VAB_matrix(3,3)]'; %L = 2
            %VAB = [xx
            %       xy
%                   xz
%                   yy
%                   yz
%                   zz];
    elseif (size(VAB_matrix,2) == 6) %For now (6/01/2017) I won't use more than f functions, so I don't think I need more for now
        %[Xxx,Xxy,Xxz,Xyy,Xyz,Xzz;
        % Yxx,Yxy,Yxz,Yyy,Yyz,Yzz;
        % Zxx,Zxy,Zxz,Zyy,Zyz,Zzz];
        VAB_L1_L2_order = [VAB_matrix(1,1:6)';VAB_matrix(2,4:6)';VAB_matrix(3,6)]'; %L = 3
%     elseif (size(VAB_matrix,2) == 10)
%         %[Xxxx,Xxxy,Xxxz,Xxyy,Xxyz,Xxzz,Xyyy,Xyyz,Xyzz,Xzzz;
%         % Yxxx,Yxxy,Yxxz,Yxyy,Yxyz,Yxzz,Yyyy,Yyyz,Yyzz,Yzzz;
%         % Zxxx,Zxxy,Zxxz,Zxyy,Zxyz,Zxzz,Zyyy,Zyyz,Zyzz,Zzzz];
%         VAB = [VAB_matrix(1,1:10)';VAB_matrix(2,7:10)';VAB_matrix(3,10)]'; %L = 4, (L+1)*(L+2)/2 = 5*6/2 = 15
        
        %The case up to L = 4 (g-functions) tells me the following:
        %VAB = [VAB_matrix(1,1:end)';VAB_matrix(2,end-L+1:end)';VAB_matrix(3,end)];
    end
    
elseif (L1 > 0 && L2 > 0)
    VAB_L1_L2_order = VAB_L1_L2_order;
end




end