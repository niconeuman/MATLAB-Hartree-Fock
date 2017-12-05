function VAB = hrr_nuc_2(basis_a,basis_b,L1,L2,Boys_Table,AL,Z,nz,unique)

RAB = -[basis_b.g(1).x0-basis_a.g(1).x0;basis_b.g(1).y0-basis_a.g(1).y0;basis_b.g(1).z0-basis_a.g(1).z0];

%L1 >= L2 always because of permutation

if L2 == 0 %[ps], [ds], [fs], etc
    VAB = contract_vrr_nuc_2(basis_a,basis_b,L1,L2,AL,Z,Boys_Table,nz);

elseif L2 == 1 %[pp], [dp], [fp], etc
    
    VABm1 = hrr_nuc_2(basis_a,basis_b,L1,L2-1,Boys_Table,AL,Z,nz,unique); %[ps],[ds],[fs], etc
    VAp1Bm1 = hrr_nuc_2(basis_a,basis_b,L1+1,L2-1,Boys_Table,AL,Z,nz,unique); %[ds],[fs],[gs], etc
    
    %Case [p|p] needs [p|s] and [d|s]
    term1 = int_expand(VABm1,RAB,L2-1,2);
    term2 = int_reshape(VAp1Bm1,L1+1,L2-1,1,unique);
        
    VAB = term1+term2;
    
elseif L2 == 2
    VABm1 = hrr_nuc_2(basis_a,basis_b,L1,L2-1,Boys_Table,AL,Z,nz,unique); %[ps],[ds],[fs], etc
    VAp1Bm1 = hrr_nuc_2(basis_a,basis_b,L1+1,L2-1,Boys_Table,AL,Z,nz,unique); %[ds],[fs],[gs], etc
    
    %Case [p|p] needs [p|s] and [d|s]
    term1 = int_expand(VABm1,RAB,L2-1,2);
    term2 = int_reshape(VAp1Bm1,L1+1,L2-1,1,unique);
        
    VAB = term1+term2;
end

end