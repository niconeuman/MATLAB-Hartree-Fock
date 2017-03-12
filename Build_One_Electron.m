function [S,T] = Build_One_Electron(basis,Shell_Doublets,NShell_Doublets)

%This function needs to move over shell doublets
%December 30th 2016
%When the basis functions involve only s and p functions, it will call
%one_electron(g1,g2)
%But when there are higher angular momentum functions, it will call
%one_electron_highL(g1,L1,g2,L2)
%this function will recursively call itself to calculate the overlap and
%kinetic energy matrix elements

Ncont = Shell_Doublets(end,2); %This is the final mu_end
S = zeros(Ncont,Ncont);
T = zeros(Ncont,Ncont);

for t = 1:NShell_Doublets
%This loop has the following characteristics. The t row of the
%Shell_Doublets matrix contains
%Shell_Doublets(t,:) = [mu_begin mu_end a nu_begin nu_end b]

%a and b run over basis functions, and are the same for all cartesian
%components of functions with angular momentum larger than zero.
%mu_begin and mu_end, etc, indicate the matrix indexes in the S and T
%matrices that will be filled by the output of the goverlap function, which
%will be generally a matrix.
mu_begin = Shell_Doublets(t,1);
mu_end = Shell_Doublets(t,2);
nu_begin = Shell_Doublets(t,4);
nu_end = Shell_Doublets(t,5);
a = Shell_Doublets(t,3);
b = Shell_Doublets(t,6);

    if (basis{a}.L < 2 && basis{b}.L < 2)
        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                [Sab,Tab] = one_electron(basis{a}.g(nba),basis{b}.g(nbb));
                
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                T(mu_begin:mu_end,nu_begin:nu_end) = T(mu_begin:mu_end,nu_begin:nu_end) + Tab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;              
            end
        end
    else
        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                [Sab,Tab] = one_electron_highL(basis{a}.g(nba),basis{a}.L,basis{b}.g(nbb),basis{b}.L);
                
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                T(mu_begin:mu_end,nu_begin:nu_end) = T(mu_begin:mu_end,nu_begin:nu_end) + Tab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;              
            end
        end
    end

end

end