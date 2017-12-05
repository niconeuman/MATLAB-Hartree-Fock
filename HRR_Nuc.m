function [VAB_L1_L2] = HRR_Nuc(RAB,V_L1_p1_L2_m1,V_L1_L2_m1,L1,L2)

Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
Dim2 = (L2+1)*(L2+2)/2;

%RAB = row vector (B-A);

if L1 == 1 && L2 == 1  %[pp]
    %Althoug the shape of the shells generated by
    %recursive_electron_nuclear_3 are column vectors, it is simpler to call
    %this function and then reorder those elements in a matrix form,
    %copying redundant elements such as yx = xy, yxx = xxy, etc
    
    %pp           =  ds (3 x 3) +    ps * RAB
    
    %[xx xy xz                
    % yx yy yz    =               
    % zx zy zz]                 
    
    V_ps = V_L1_L2_m1';       
    V_ds = [V_L1_p1_L2_m1(1:3),V_L1_p1_L2_m1([2 4 5]),V_L1_p1_L2_m1([3 5 6])];
    
    VAB_L1_L2 = diag(diag(V_ds)) + diag(V_ps.*RAB);
    
    
    
elseif L1 == 2 && L2 == 1  %[dp]
    %Althoug the shape of the shells generated by
    %recursive_electron_nuclear_3 are column vectors, it is simpler to call
    %this function and then reorder those elements in a matrix form,
    %copying redundant elements such as yx = xy, yxx = xxy, etc
    
    %dp (6 x 3)   =  fs (6 x 3) +    ds * RAB
    
    %fs = [xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz]' -> 
    %[fs(1:6) fs([2 4 5 7 8 9]) fs([3 5 6 8 9 10])]
   
    
    %[xxX xxY xxZ                
    % xyX                        
    % xzX         =          
    % yyX
    % yzX
    % zzX    zzZ]
    
    V_ds = V_L1_L2_m1;
    V_fs = [V_L1_p1_L2_m1(1:6),V_L1_p1_L2_m1([2 4 5 7 8 9]),V_L1_p1_L2_m1([3 5 6 8 9 10])];
    %Note, the indexes in the third column are the indices in the second +
    %1
    
    VAB_L1_L2 = V_fs + V_ds*RAB;
elseif L1 == 2 && L2 == 2 %[dd]
    
    %dd (6 x 6)   =  fp (6 x 6) +    dp (6 x 3) * RAB
    
    %fp = [xxxX xxxY xxxZ
%          xxyX xxyY xxyZ
%          xxzX xxzY xxzY
%          xyyX xyyY
%          xyzX xyzY
%          xzzX xzzY
%          yyyX yyyY
%          yyzX yyzY
%          yzzX yzzY
%          zzzX zzzY      ] -> 
    %[fp(1:6,:) fp([2 4 5 7 8 9],2:3) fp([3 5 6 8 9 10],3)]
%   V_fp  = xx xy xz yy yz zz
%     xx
%     xy
%     xz
%     yy
%     yz
%     zz
    
    %[xxX xxY xxZ                
    % xyX                        
    % xzX         =          
    % yyX
    % yzX
    % zzX    zzZ]
    V_dp = V_L1_L2_m1;
    V_fp = [V_L1_p1_L2_m1(1:6,:),V_L1_p1_L2_m1([2 4 5 7 8 9],2:3),V_L1_p1_L2_m1([3 5 6 8 9 10],3)];
    
    VAB_L1_L2 = V_fp + [RAB(1)*V_dp,RAB(2)*V_dp(:,2:3);RAB(3)*V_dp(:,3)];
    
end
end