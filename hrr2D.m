function gabcd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique)

RAB = -[basis_b.g(1).x0-basis_a.g(1).x0;basis_b.g(1).y0-basis_a.g(1).y0;basis_b.g(1).z0-basis_a.g(1).z0];
%RCD = -[basis_d.g(1).x0-basis_c.g(1).x0;basis_d.g(1).y0-basis_c.g(1).y0;basis_d.g(1).z0-basis_c.g(1).z0];

if L2 == 0 %[ps|ss], [ds|ss], [fs|ss], etc
    if (L3 == 0 && L4 == 0)
        gabcd = contract_vrr(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz);
    elseif (L3 > 0 && L4 == 0)
        if (L3 > L1)
           gabcd = contract_ascs(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
        else
           gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
        end
        
        %gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
        %gabcd = permute(gabcd, [1 3 2 4]);
    elseif (L3 == 0 && L4 > 0)
        if (L4 > L1)
           gabcd = contract_ascs(basis_d,basis_c,basis_a,basis_b,L4,L3,L1,L2,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
        else
           gabcd = contract_ascs(basis_a,basis_b,basis_d,basis_c,L1,L2,L4,L3,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
        end
        
        %gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
        
        gabcd = permute(gabcd, [1 2 4 3]);
    elseif (L3 > 0 && L4 > 0)
        %I'm calculating a [LaS|LcLd] type integral, so I need to permute
        %the bra and ket, do hrr2D, and permute back
        if (L4 > L3)
           gabcd = hrr2D(basis_d,basis_c,basis_a,basis_b,L4,L3,L1,L2,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [2 1 3 4]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
        else
           gabcd = hrr2D(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
        end
        
        
        gabcd = permute(gabcd, [3 4 1 2]);
    end
elseif L2 == 1 %[pp|ss], [dp|ss], [fp|ss], etc
    
    gabm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1,L2-1,L3,L4,Boys_Table,nz,unique); %[ps],[ds],[fs], etc
    gap1bm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1+1,L2-1,L3,L4,Boys_Table,nz,unique); %[ds],[fs],[gs], etc
    
    %Case [pp|ss] needs [ps|ss] and [ds|ss]
    term1 = int_expand(gabm1cd,RAB,L2-1,2);
    term2 = int_reshape(gap1bm1cd,L1+1,L2-1,1,unique);
        
    gabcd = term1+term2;
    
elseif L2 == 2
    
    gabm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1,L2-1,L3,L4,Boys_Table,nz,unique); %[ps],[ds],[fs], etc
    gap1bm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1+1,L2-1,L3,L4,Boys_Table,nz,unique); %[ds],[fs],[gs], etc
    
    %Case [pp|ss] needs [ps|ss] and [ds|ss]
    term1 = int_expand(gabm1cd,RAB,L2-1,2);
    term2 = int_reshape(gap1bm1cd,L1+1,L2-1,1,unique);
        
    gabcd = term1+term2;
    
end





end