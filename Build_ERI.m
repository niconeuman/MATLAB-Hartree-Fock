function gabcd = Build_ERI(basis,Shells,NShells,Boys_Table)

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
gabcd = cell(Ncont,Ncont);


for t = 1:NShells
    mu_begin = Shells(t,1);
    mu_end = Shells(t,2);
    nu_begin = Shells(t,4);
    nu_end = Shells(t,5);
    ka_begin = Shells(t,7);
    ka_begin = Shells(t,8);
    la_begin = Shells(t,10);
    la_begin = Shells(t,11);
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
    
    if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals
        for na = 1:basis_a.n
    
        end
    elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
    elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
    elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals    
    elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
    elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals
    elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals
    elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals
    elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals
    elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals
    elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals    
    elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals
    elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals
    elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals
    elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals    
    elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals
    end
end











end