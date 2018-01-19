function gabcd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique)

RAB = -[basis_b.g(1).x0-basis_a.g(1).x0;basis_b.g(1).y0-basis_a.g(1).y0;basis_b.g(1).z0-basis_a.g(1).z0];
RCD = -[basis_d.g(1).x0-basis_c.g(1).x0;basis_d.g(1).y0-basis_c.g(1).y0;basis_d.g(1).z0-basis_c.g(1).z0];

if L2 == 0 %[ps|ss], [ds|ss], [fs|ss], etc
    if (L1 == 0)
        if (L3 > 0 && L4 == 0)
            [gabcd,~] = contract_vrr_2(basis_c,basis_d,basis_a,basis_b,L3,0,0,0,Boys_Table);
            %gabcd = hrr2D(basis_b,basis_a,basis_c,basis_d,L2,L1,L3,L4,Boys_Table,nz,unique);
            gabcd = permute(gabcd, [3 4 1 2]);         
        elseif (L3 == 0 && L4 > 0)
            [gabcd,~] = contract_vrr_2(basis_d,basis_c,basis_a,basis_b,L4,0,0,0,Boys_Table);
            %gabcd = hrr2D(basis_b,basis_a,basis_c,basis_d,L2,L1,L3,L4,Boys_Table,nz,unique);
            gabcd = permute(gabcd, [3 4 1 2]);
            gabcd = permute(gabcd, [1 2 4 3]);
        elseif (L3 > 0 && L4 > 0)
                if (L4 > L3)
                   gabcd = hrr2D(basis_d,basis_c,basis_a,basis_b,L4,L3,L1,L2,Boys_Table,nz,unique);
                   gabcd = permute(gabcd, [2 1 3 4]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
                else
                   gabcd = hrr2D(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
                end
                gabcd = permute(gabcd, [3 4 1 2]);
        end
    elseif (L1 > 0)
        if L3 == 0 && L4 == 0
                [gabcd,~] = contract_vrr_2(basis_a,basis_b,basis_c,basis_d,L1,0,0,0,Boys_Table);
        elseif (L3 > 0 && L4 == 0)
                if (L3 > L1)
                   gabcd = contract_ascs_2(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
                   gabcd = permute(gabcd, [1 3 2 4]);
                   gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
                else
                   gabcd = contract_ascs_2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
                   gabcd = permute(gabcd, [1 3 2 4]);
                end
        elseif (L3 == 0 && L4 > 0)
                if (L4 > L1)
                   gabcd = contract_ascs_2(basis_d,basis_c,basis_a,basis_b,L4,L3,L1,L2,Boys_Table,nz,unique);
                   gabcd = permute(gabcd, [1 3 2 4]);
                   gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
                else
                   gabcd = contract_ascs_2(basis_a,basis_b,basis_d,basis_c,L1,L2,L4,L3,Boys_Table,nz,unique);
                   gabcd = permute(gabcd, [1 3 2 4]);
                end
                gabcd = permute(gabcd, [1 2 4 3]);
        else %L3 > 0 && L4 > 0
                %gascs,gap1scm1s,gam1scm1s,gascm1s,gascm2s
                %In the future I might have to take care of the case where
                %d.L > c.L
                [permutedg2010,~,~,~,~] = contract_ascs_2(basis_c,basis_d,basis_b,basis_a,2,0,1,0,Boys_Table,nz,unique);
                [permutedg1010,~,~,~,~] = contract_ascs_2(basis_c,basis_d,basis_b,basis_a,1,0,1,0,Boys_Table,nz,unique);
                permutedg2010 = permute(permutedg2010,[1 3 2 4]);
                permutedg1010 = permute(permutedg1010,[1 3 2 4]);

                gascm1s = permute(permutedg1010,[1 2 4 3]); %To obtain cdab
                gap1scm1s = permute(permutedg2010,[1 2 4 3]); %To obtain cdab
                
                disp('This case occurred');
                permuted_term1 = int_expand(gascm1s,RCD,L2-1,2); %RCD, not RAB
                permuted_term2 = int_reshape(gap1scm1s,L1+1,L2-1,1,unique);
        
                permuted_gabcd = permuted_term1+permuted_term2;
                gabcd = permute(permuted_gabcd,[3 4 1 2]);
        end
%             if (L1 > 0 && L3 == 0 && L4 == 0)
%                 gabcd = contract_vrr(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz);
%             elseif (L1 == 0 && L3 == 0 && L4 == 0)
%                 %Replace this for a newer version of contract_ssss
%                 %gabcd = contract_ssss(basis_a,basis_b,basis_c,basis_d,Boys_Table);
%                 disp('ssss integrals should not occur in hrr2D');
%             elseif (L3 > 0 && L4 == 0)
%                 if (L3 > L1)
%                    gabcd = contract_ascs_2(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
%                    gabcd = permute(gabcd, [1 3 2 4]);
%                    gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
%                 else
%                    gabcd = contract_ascs_2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
%                    gabcd = permute(gabcd, [1 3 2 4]);
%                 end
% 
%                 %gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
%                 %gabcd = permute(gabcd, [1 3 2 4]);
%             elseif (L3 == 0 && L4 > 0)
%                 if (L4 > L1)
%                    gabcd = contract_ascs_2(basis_d,basis_c,basis_a,basis_b,L4,L3,L1,L2,Boys_Table,nz,unique);
%                    gabcd = permute(gabcd, [1 3 2 4]);
%                    gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
%                 else
%                    gabcd = contract_ascs_2(basis_a,basis_b,basis_d,basis_c,L1,L2,L4,L3,Boys_Table,nz,unique);
%                    gabcd = permute(gabcd, [1 3 2 4]);
%                 end
% 
%                 %gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique);
% 
%                 gabcd = permute(gabcd, [1 2 4 3]);
%             elseif (L3 > 0 && L4 > 0)
%                 %I'm calculating a [LaS|LcLd] type integral, so I need to permute
%                 %the bra and ket, do hrr2D, and permute back
%                 if (L4 > L3)
%                    gabcd = hrr2D(basis_d,basis_c,basis_a,basis_b,L4,L3,L1,L2,Boys_Table,nz,unique);
%                    gabcd = permute(gabcd, [2 1 3 4]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
%                 else
%                    gabcd = hrr2D(basis_c,basis_d,basis_a,basis_b,L3,L4,L1,L2,Boys_Table,nz,unique);
%                 end
% 
% 
%                 gabcd = permute(gabcd, [3 4 1 2]);
%             end
    end
elseif L2 == 1 %[pp|ss], [dp|ss], [fp|ss], etc

    if (L1 == 0 && L3 == 0 && L4 == 0)
            [gabcd,~] = contract_vrr_2(basis_b,basis_a,basis_c,basis_d,L2,L1,L3,L4,Boys_Table);
            %gabcd = hrr2D(basis_b,basis_a,basis_c,basis_d,L2,L1,L3,L4,Boys_Table,nz,unique);
            gabcd = permute(gabcd, [2 1 3 4]);
    elseif (L1 == 0 && L3 > 0 && L4 == 0)
        if (L3 > L2)
           gabcd = contract_ascs_2(basis_c,basis_d,basis_b,basis_a,L3,L4,L2,L1,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
           gabcd = permute(gabcd, [2 1 3 4]);
        else
           gabcd = contract_ascs_2(basis_b,basis_a,basis_c,basis_d,L2,L1,L3,L4,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [2 1 3 4]);
        end
    elseif (L1 == 0 && L3 == 0 && L4 > 0)
        if (L4 > L2)
           gabcd = contract_ascs_2(basis_d,basis_c,basis_b,basis_a,L4,L3,L2,L1,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
           gabcd = permute(gabcd, [2 1 4 3]);
        else
           gabcd = contract_ascs_2(basis_b,basis_a,basis_d,basis_c,L2,L1,L3,L4,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [1 3 2 4]);
           gabcd = permute(gabcd, [2 1 4 3]);
        end
    elseif (L1 == 0 && L3 > 0 && L4 > 0)
        if (L4 > L3)
           gabcd = hrr2D(basis_d,basis_c,basis_b,basis_a,L4,L3,L2,L1,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [3 4 1 2]); %remember ascs produces a 2D matrix. So permutation of L3 and L1 transposes the matrix
           gabcd = permute(gabcd, [2 1 4 3]);
        else
           gabcd = hrr2D(basis_c,basis_d,basis_b,basis_a,L3,L4,L2,L1,Boys_Table,nz,unique);
           gabcd = permute(gabcd, [3 4 1 2]);
           gabcd = permute(gabcd, [2 1 3 4]);
        end        
    else %(L1 ==1 && L3,L4 not both 1, because the [pp|pp] case is handled by hrr4D)
            %gabm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1,L2-1,L3,L4,Boys_Table,nz,unique); %[ps],[ds],[fs], etc
            %gap1bm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1+1,L2-1,L3,L4,Boys_Table,nz,unique); %[ds],[fs],[gs], etc
            if L3 == 0 && L4 == 0
                [g2000,g1000] = contract_vrr_2(basis_a,basis_b,basis_c,basis_d,2,0,0,0,Boys_Table);
                gascm1s = g1000;
                gap1scm1s = g2000;
            term1 = int_expand(gascm1s,RAB,L2-1,2);
            term2 = int_reshape(gap1scm1s,L1+1,L2-1,1,unique);
        
            gabcd = term1+term2;
            elseif L3 == 1 && L4 == 0
                %gascs,gap1scm1s,gam1scm1s,gascm1s,gascm2s
                [g2010,g3000,g1000,g2000,~] = contract_ascs_2(basis_a,basis_b,basis_c,basis_d,2,0,1,0,Boys_Table,nz,unique);
                [g1010,g2000,g0000,g1000,~] = contract_ascs_2(basis_a,basis_b,basis_c,basis_d,1,0,1,0,Boys_Table,nz,unique);
                g2010 = permute(g2010,[1 3 2 4]);
                g1010 = permute(g1010,[1 3 2 4]);
                gascm1s = g1010;
                gap1scm1s = g2010;
            term1 = int_expand(gascm1s,RAB,L2-1,2);
            term2 = int_reshape(gap1scm1s,L1+1,L2-1,1,unique);
        
            gabcd = term1+term2;    
            elseif L3 == 0 && L4 == 1
                %gascs,gap1scm1s,gam1scm1s,gascm1s,gascm2s
                [g2010,g3000,g1000,g2000,~] = contract_ascs_2(basis_a,basis_b,basis_d,basis_c,2,0,1,0,Boys_Table,nz,unique);
                [g1010,g2000,g0000,g1000,~] = contract_ascs_2(basis_a,basis_b,basis_d,basis_c,1,0,1,0,Boys_Table,nz,unique);
                g2010 = permute(g2010,[1 3 2 4]);
                g1010 = permute(g1010,[1 3 2 4]);
                gascm1s = permute(g1010,[1 2 4 3]);
                gap1scm1s = permute(g2010,[1 2 4 3]);
            term1 = int_expand(gascm1s,RAB,L2-1,2);
            term2 = int_reshape(gap1scm1s,L1+1,L2-1,1,unique);
        
            gabcd = term1+term2;    
            end    
            %Case [pp|ss] needs [ps|ss] and [ds|ss]

     
    end

    
elseif L2 == 2
    disp('This should not occur for [pp|pp] integrals');
    gabm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1,L2-1,L3,L4,Boys_Table,nz,unique); %[ps],[ds],[fs], etc
    gap1bm1cd = hrr2D(basis_a,basis_b,basis_c,basis_d,L1+1,L2-1,L3,L4,Boys_Table,nz,unique); %[ds],[fs],[gs], etc
    
    %Case [pp|ss] needs [ps|ss] and [ds|ss]
    term1 = int_expand(gabm1cd,RAB,L2-1,2);
    term2 = int_reshape(gap1bm1cd,L1+1,L2-1,1,unique);
        
    gabcd = term1+term2;
    
end





end