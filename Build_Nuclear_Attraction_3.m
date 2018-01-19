function [VAB,nAtom] = Build_Nuclear_Attraction_3(basis,Shell_List,AL,Z,Boys_Table)
nAtom = size(Z,2);
Ncont = Shell_List(end,2);
nb = size(basis,1);
VAB = zeros(Ncont,Ncont);



nz = cell(1,10);
nz{2} = [1;0;0;1;0;1];
nz{3} = [2;1;1;0;0;0;2;1;0;2];
nz{4} = [3;2;2;1;1;1;0;0;0;0;3;2;1;0;3];
nz{5} = [4;3;3;2;2;2;1;1;1;1;0;0;0;0;0;4;3;2;1;0;4];
nz{6} = [5;4;4;3;3;3;2;2;2;2;1;1;1;1;1;0;0;0;0;0;0;5;4;3;2;1;0;5];
nz{7} = [6;5;5;4;4;4;3;3;3;3;2;2;2;2;2;1;1;1;1;1;1;0;0;0;0;0;0;0;6;5;4;3;2;1;0;6];
nz{8} = [7;6;6;5;5;5;4;4;4;4;3;3;3;3;3;2;2;2;2;2;2;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;7;6;5;4;3;2;1;0;7];
nz{9} = [8;7;7;6;6;6;5;5;5;5;4;4;4;4;4;3;3;3;3;3;3;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;8;7;6;5;4;3;2;1;0;8];
nz{10} = [9;8;8;7;7;7;6;6;6;6;5;5;5;5;5;4;4;4;4;4;4;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;9;8;7;6;5;4;3;2;1;0;9];

unique = [2, 4,5, 7,8,9, 11,12,13,14, 16,17,18,19,20, 22,23,24,25,26,27,  29,30,31,32,33,34,35, 37,38,39,40,41,42,43,44, 46,47,48,49,50,51,52,53,54, 56,57,58,59,60,61,62,63,64,65 ]; %up to [ms|ss]

for a = 1:nb
    mu_begin = Shell_List(a,1);
    mu_end = Shell_List(a,2);
    %aa = Shells_List(a,2);
    basis_a = basis{a};
    for b = 1:a
        nu_begin = Shell_List(b,1);
        nu_end = Shell_List(b,2);
        %ab = Shells_List(b,2);
        basis_b = basis{b};
        
        
        if (basis_a.L == 0 && basis_b.L == 0) %[S|eN|S] integrals
            VAB(mu_begin:mu_end,nu_begin:nu_end) = contract_ss_Nuc(basis_a,basis_b,Boys_Table,AL,Z);
            VAB(nu_begin:nu_end,mu_begin:mu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end);
        elseif (basis_a.L > 0 && basis_b.L == 0) %[L|eN|S] integrals
            VAB(mu_begin:mu_end,nu_begin:nu_end) = contract_vrr_nuc_2(basis_a,basis_b,basis_a.L,basis_b.L,AL,Z,Boys_Table,nz);
            VAB(nu_begin:nu_end,mu_begin:mu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end);
        elseif (basis_a.L == 0 && basis_b.L > 0) %[S|eN|L] integrals
            permuted_VAB = contract_vrr_nuc_2(basis_b,basis_a,basis_b.L,basis_a.L,AL,Z,Boys_Table,nz);
            VAB(mu_begin:mu_end,nu_begin:nu_end) = permuted_VAB';
            VAB(nu_begin:nu_end,mu_begin:mu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end);
        else %(basis_a.L > 0 && basis_b.L > 0)
            if (basis_b.L > basis_a.L)
                VABtemp = hrr_nuc_2(basis_b,basis_a,basis_b.L,basis_a.L,Boys_Table,AL,Z,nz,unique);
                VABtemp = VABtemp';
            else
                VABtemp = hrr_nuc_2(basis_a,basis_b,basis_a.L,basis_b.L,Boys_Table,AL,Z,nz,unique);
            end
            VAB(mu_begin:mu_end,nu_begin:nu_end) = VABtemp;
            VAB(nu_begin:nu_end,mu_begin:mu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end);
        end
        
    end
end

