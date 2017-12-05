function gabcd = Build_ERI(basis,Shells,NShells,Boys_Table,pair_data)

%February 13th 2017
%This function calculates the two-electron repulsion matrix for
%conventional Hartree-Fock, from a list of Shell quartets. Eventually
%separate programs could be used for SSSS type integrals (listed in
%SSSS_Shells), (LaLb|LcLd) shells with some degree of contraction (listed
%in Contracted_HighL_Shells) and uncontracted (LaLb|LcLd) (listed in
%Uncontracted_HighL_Shells).

Ncont = Shells(end,2); %This is the final mu_end
%In Matlab there are no 3D or 4D sparse arrays (yet).
%What seems to be easiest is to make a 2D cell array (which is by
%definition sparse) and fill it with 2D sparse arrays. This second step
%maybe should be performed in the loop over a and b (just to save a 2D
%loop).
%A 100^4 matrix of zeros takes 800 Mb
%A 100x100 cell array of 100x100 sparse matrices takes 9.36 Mb
gabcd = zeros(Ncont,Ncont,Ncont,Ncont);

%This will give me the nonzero elements that contribute in the [a-2s|ss]
%terms in the vertical recursion relations.
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

for t = 1:NShells
    mu_begin = Shells(t,1);
    mu_end = Shells(t,2);
    nu_begin = Shells(t,4);
    nu_end = Shells(t,5);
    ka_begin = Shells(t,7);
    ka_end = Shells(t,8);
    la_begin = Shells(t,10);
    la_end = Shells(t,11);
    a = Shells(t,3);
    b = Shells(t,6);
    c = Shells(t,9);
    d = Shells(t,12);
    
    basis_a = basis{a};
    basis_b = basis{b};
    basis_c = basis{c};
    basis_d = basis{d};
    
    %Now I need to separate into the different cases
    %There are 16 types.
    ab_data = pair_data{a,b};
    cd_data = pair_data{c,d};
    
    if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals
        gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_2(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
    elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
         gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_vrr(basis_a,basis_b,basis_c,basis_d,Boys_Table,nz);
     
    elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
         permuted_LSSS = contract_vrr(basis_b,basis_a,basis_c,basis_d,Boys_Table,nz); %b and a are permuted
         gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [2 1 3 4]);
    elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
        permuted_LSSS = contract_vrr(basis_c,basis_b,basis_a,basis_d,Boys_Table,nz); %b and a are permuted
        gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [3 2 1 4]);
    elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
        permuted_LSSS = contract_vrr(basis_d,basis_b,basis_c,basis_a,Boys_Table,nz); %b and a are permuted
        gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [4 2 3 1]);
    
         
    elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals
        LaSLcS = contract_ascs(basis_a,basis_b,basis_c,basis_d,Boys_Table,nz,unique);
        %LaSLcS is structured as a [DimA, DimC,1,1] matrix, for convenience in matrix operations.
        LaSLcS = permute(LaSLcS, [1 3 2 4]);
        gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;
%     elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals
%     elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals
%     elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals
%     elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals
%     elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals    
%     elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals
%     elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals
%     elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals
%     elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals    
%     elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals

    else
        gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);
%
    end
end











end