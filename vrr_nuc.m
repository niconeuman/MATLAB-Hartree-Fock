function [Vab0,Vab1] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1,L2,order,Boys_Table,nz)

Dim1 = (L1+1)*(L1+2)/2;

if L1 == 0
q = a*b/p;
Kab = exp(-q*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
prefactor = -1*Z*(2*pi/p)*Kab;
Vab0 = prefactor*Interpolated_Boys(order,p*RPC2,Boys_Table);
Vab1 = prefactor*Interpolated_Boys(order+1,p*RPC2,Boys_Table);
    


elseif L1 == 1
    [Vam1b0,~] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order,Boys_Table,nz);
    [Vam1b1,Vam1b2] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order+1,Boys_Table,nz);
    
    term1_0 = int_expand(Vam1b0,RPA,L1-1,1);
    term2_0 = int_expand(Vam1b1,RPC,L1-1,1);
    term1_1 = int_expand(Vam1b1,RPA,L1-1,1);
    term2_1 = int_expand(Vam1b2,RPC,L1-1,1);
    
    Vab0 = term1_0-term2_0;
    Vab1 = term1_1-term2_1;
    
elseif L1 == 2
    [Vam2b0,~] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,order,Boys_Table,nz);
    [Vam2b1,Vam2b2] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,order+1,Boys_Table,nz);
    [Vam1b0,~] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order,Boys_Table,nz);
    [Vam1b1,Vam1b2] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order+1,Boys_Table,nz);

    term1_0 = int_expand(Vam1b0,RPA,L1-1,1);
    term2_0 = int_expand(Vam1b1,RPC,L1-1,1);
    term3_0 = (1/2/p)*nz{L1}.*(Vam2b0-Vam2b1);
       
    
    term1_1 = int_expand(Vam1b1,RPA,L1-1,1);
    term2_1 = int_expand(Vam1b2,RPC,L1-1,1);
    term3_1 = (1/2/p)*nz{L1}.*(Vam2b1-Vam2b2);

    
    Vab0 = term1_0-term2_0+term3_0;
    Vab1 = term1_1-term2_1+term3_1;

else %(L1 > 2)
    [Vam2b0,~] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,order,Boys_Table,nz);
    [Vam2b1,Vam2b2] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-2,L2,order+1,Boys_Table,nz);
    [Vam1b0,~] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order,Boys_Table,nz);
    [Vam1b1,Vam1b2] = vrr_nuc(a,b,RAB,RPA,RPB,RPC,p,RPC2,Z,L1-1,L2,order+1,Boys_Table,nz);

    term1_0 = int_expand(Vam1b0,RPA,L1-1,1);
    term2_0 = int_expand(Vam1b1,RPC,L1-1,1);
    %Difference with case L1 == 2
    prevterm3_0 = Vam2b0-Vam2b1;
    expterm3_0 = int_expand(prevterm3_0,ones(3,1),L1-2,1);
    term3_0 = (1/2/p)*nz{L1}.*(int_expand(expterm3_0,ones(3,1),L1-1,1));
       
    
    term1_1 = int_expand(Vam1b1,RPA,L1-1,1);
    term2_1 = int_expand(Vam1b2,RPC,L1-1,1);
    %Difference with case L1 == 2
    prevterm3_1 = Vam2b1-Vam2b2;
    expterm3_1 = int_expand(prevterm3_1,ones(3,1),L1-2,1);
    term3_1 = (1/2/p)*nz{L1}.*(int_expand(expterm3_1,ones(3,1),L1-1,1));

    
    Vab0 = term1_0-term2_0+term3_0;
    Vab1 = term1_1-term2_1+term3_1;

end



end