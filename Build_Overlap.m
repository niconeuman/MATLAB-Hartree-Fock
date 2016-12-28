function [S,T] = Build_Overlap(basis,Shell_Doublets,NShell_Doublets)
%This function needs to move over shell doublets
%nb = size(basis,2);
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

        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + goverlap(basis{a}.g(nba),basis{b}.g(nbb))*basis{a}.c(nba)*basis{b}.c(nbb);
                %T(n,m) = T(n,m) + gkinetic(basis{n}.g(nba),basis{m}.g(nbb))*basis{n}.c(nba)*basis{m}.c(nbb);
            end
        end


end




end