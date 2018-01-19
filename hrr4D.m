function gabcd = hrr4D(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique)

%December 27th 2017
%This function does horizontal recursion on the bra or ket, in a predefined
%order.

if ((L1 == 0 || L2 == 0) && (L3 == 0 || L4 == 0)) %Just one of each has non-zero angular momentum
    disp('This should not be covered by hrr4D');
    if (L1 > 0 && L3 > 0)
        if (L3 > L1)
           gabcd = contract_ascs_2(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
        else
           gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
        end
        
    end
elseif (L1 == 1 && L2 == 1 && L3 == 1 && L4 == 1)
    
    [gc2d0a2b0,~,gc1d0a1b0,gc2d0a1b0,~] = contract_ascs_2(basis_c,basis_d,basis_a,basis_b,2,0,2,0,Boys_Table,nz,unique);
    %Calculate g1010,g2020,g2010,g1020
    %g2020permuted = contract_ascs_2(basis_c,basis_d,basis_a,basis_b,2,0,2,0,Boys_Table,nz,unique);
    gc2d0a2b0 = permute(gc2d0a2b0, [1 3 2 4]);
    gc1d0a1b0 = permute(gc1d0a1b0, [1 3 2 4]);
    gc2d0a1b0 = permute(gc2d0a1b0, [1 3 2 4]);
    
    ga2b0c1d0 = contract_ascs_2(basis_a,basis_b,basis_c,basis_d,2,0,1,0,Boys_Table,nz,unique);
    ga2b0c1d0 = permute(ga2b0c1d0, [1 3 2 4]);
    %g2010permuted = permute(ga2b0c1d0, [3 4 1 2]);
    gc1d0a2b0 = permute(ga2b0c1d0, [3 4 1 2]);
    
    
RAB = [basis_a.g(1).x0-basis_b.g(1).x0;basis_a.g(1).y0-basis_b.g(1).y0;basis_a.g(1).z0-basis_b.g(1).z0];
RCD = [basis_c.g(1).x0-basis_d.g(1).x0;basis_c.g(1).y0-basis_d.g(1).y0;basis_c.g(1).z0-basis_d.g(1).z0];    
    
    gc1d1a2b0_term1 = int_expand(gc1d0a2b0,RCD,0,2);
    gc1d1a2b0_term2 = int_reshape(gc2d0a2b0,2,0,1,unique);
    
    %This is still permutted
    gc1d1a2b0 = gc1d1a2b0_term1+gc1d1a2b0_term2;
    
    gc1d1a1b0_term1 = int_expand(gc1d0a1b0,RCD,0,2);
    gc1d1a1b0_term2 = int_reshape(gc2d0a1b0,2,0,1,unique);
    
    %This is still permutted
    gc1d1a1b0 = gc1d1a1b0_term1+gc1d1a1b0_term2;
    
    ga2b0c1d1 = permute(gc1d1a2b0, [3 4 1 2]);
    ga1b0c1d1 = permute(gc1d1a1b0, [3 4 1 2]);
    
    g1111term1 = int_expand(ga1b0c1d1,RAB,0,2);
    g1111term2 = int_reshape(ga2b0c1d1,2,0,1,unique);
    
    gabcd = g1111term1 + g1111term2;
end





end