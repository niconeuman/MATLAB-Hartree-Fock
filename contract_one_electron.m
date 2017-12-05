function [S,T,Pcell,pcell] = contract_one_electron(basis_a,basis_b,nz,unique)

%requires basis_a and basis_b

L1 = basis_a.L;
L2 = basis_b.L;

Dim1 = (L1+1)*(L1+2)/2;
Dim2 = (L2+1)*(L2+2)/2;
    
tempP = [0;0;0];
tempp = 0;

S = zeros(Dim1,Dim2);
T = zeros(Dim1,Dim2);
        for nba = 1:basis_a.n
            for nbb = 1:basis_b.n
                [Sab,Tab,Pab,p] = one_electron_2(basis_a.g(nba),basis_b.g(nbb));
                %[Sab,Tab] = one_electron(basis{a}.g(nba),basis{b}.g(nbb));
                
                S(Dim1,Dim2) = S(Dim1,Dim2) + Sab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                T(Dim1,Dim2) = T(Dim1,Dim2) + Tab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                tempP = tempP + Pab*Sab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N; %maybe later I need to add contraction and normalization coefficient 
                tempp = tempp + p*Sab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
            end
        end
        %For matrix elements larger than (1x1) I need to modify this in
        %some way.
        %July 11th 2017
        %P is always a 3x1 vector, and p is always a scalar. The problem is
        %that the overlap matrix will change size depending on L1 and L2,
        %but the whole shell doublet will have the same P and p. So what
        %probably should be done is to output an S-type overlap matrix
        %element (from which the whole shell doublet was obtained), and use
        %this to normalize P and p.
        
        Pcell = tempP/S(Dim1,Dim2);
        pcell = tempp/S(Dim1,Dim2);



end