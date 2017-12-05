function [VAB_L1_L2_order,VAB_L1_L2_order_p1,bufferInts] = VRR_for_HRR_Nuc(a,b,RPA,RPB,RPC,p,RPC2,Z,L1,L2,Lmax,order,Boys_Table,bufferInts)

%January 31st 2017
%This function calculates nuclear attraction matrices of the type [L1s] for
%use with 

%Data I need
%L1 = total angular momentum of the first basis function
%L2 = total angular momentum of the second basis function

%According to Helgaker's book, p. 378
Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
Dim2 = (L2+1)*(L2+2)/2;

VAB_L1_L2_order = zeros(Dim1,Dim2);
VAB_L1_L2_order_p1 = zeros(Dim1,Dim2);


if (L1 == 0 && L2 == 0)
    %In each block I need to make the function check if the bufferInts{ind}
    %contains a double. If it does, that will be the output of the
    %function. If it doesn't, it will calculate it.
    
    Kab = gprod(RPA(1),RPA(2),RPA(3),a,RPB(1),RPB(2),RPB(3),b);
    
    %31 jan 2017. Consider taking out the -1*Z*(2*pi) out to the outside
    %loop in Build_Nuclear_Attraction.
    
    for order = 0:Lmax
    Vss_m = -1*Z*(2*pi/p)*Interpolated_Boys(order,p*RPC2,Boys_Table)*Kab;
    %Generally (ind=L1*L2max*(Lmax)+L2*(Lmax)+(m+1))
    ind = order+1; %As L1 and L2 are zero, I only need the final part of the compound index
    bufferInts{ind} = Vss_m;
    end
    VAB_L1_L2_order = bufferInts{1};
    VAB_L1_L2_order_p1 = bufferInts{2};
    
elseif (L1 > 0 && L2 == 0)
    %This line of integrals will always be a COLUMN vector
    [VAB_L1_m1_order,VAB_L1_m1_order_p1] = recursive_electron_nuclear_3(a,b,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,Lmax,order,Boys_Table,bufferInts);
    
    switch L1-2
        case -1
            VAB_L1_L2_order = [RPA(1)*VAB_L1_m1_order;RPA(2)*VAB_L1_m1_order(end-L1+1:end);RPA(3)*VAB_L1_m1_order(end)]+...
                        +[RPC(1)*VAB_L1_m1_order_p1;RPC(2)*VAB_L1_m1_order_p1(end-L1+1:end);RPC(3)*VAB_L1_m1_order_p1(end)];
        case 0
            [VAB_L1_m2_order,VAB_L1_m2_order_p1] = recursive_electron_nuclear_3(a,b,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,Lmax,order,Boys_Table,bufferInts);
        %The outputs are [ss](m), and the final integral VAB_L1_L2_order =
        %[ds]
            VAB_L1_m2_order = VAB_L1_m2_order*ones(3,1);
            VAB_L1_m2_order_p1 = VAB_L1_m2_order_p1*ones(3,1);
            VAB_L1_L2_order = [RPA(1)*VAB_L1_m1_order;RPA(2)*VAB_L1_m1_order(end-L1+1:end);RPA(3)*VAB_L1_m1_order(end)]+...
                        +[RPC(1)*VAB_L1_m1_order_p1;RPC(2)*VAB_L1_m1_order_p1(end-L1+1:end);RPC(3)*VAB_L1_m1_order_p1(end)]...
                        +1/2/p*(L1-1)*([VAB_L1_m2_order-VAB_L1_m2_order_p1;VAB_L1_m2_order(end-L1+1:end)-VAB_L1_m2_order_p1(end-L1+1:end);VAB_L1_m2_order(end)-VAB_L1_m2_order_p1(end)]);

        otherwise
            [VAB_L1_m2_order,VAB_L1_m2_order_p1] = recursive_electron_nuclear_3(a,b,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,Lmax,order,Boys_Table,bufferInts);
        %The outputs are [ps](m), and the final integral is [fs] or
        %[ds](m) and the final integrals [gs], etc
        %With this I generate an array the same size as VAB_L1_m1_order(_p1)
            VAB_L1_m2_order = [VAB_L1_m2_order;VAB_L1_m2_order(end-(L1-1)+1:end);VAB_L1_m2_order(end)];
            
            VAB_L1_L2_order = [RPA(1)*VAB_L1_m1_order;RPA(2)*VAB_L1_m1_order(end-L1+1:end);RPA(3)*VAB_L1_m1_order(end)]+...
                        +[RPC(1)*VAB_L1_m1_order_p1;RPC(2)*VAB_L1_m1_order_p1(end-L1+1:end);RPC(3)*VAB_L1_m1_order_p1(end)]...
                        +1/2/p*(L1-1)*([VAB_L1_m2_order-VAB_L1_m2_order_p1;VAB_L1_m2_order(end-L1+1:end)-VAB_L1_m2_order_p1(end-L1+1:end);VAB_L1_m2_order(end)-VAB_L1_m2_order_p1(end)]);
    end
            

                    

end



end