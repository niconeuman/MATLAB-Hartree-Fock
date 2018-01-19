function [S,T,Pcell,pcell,S00] = contract_one_electron_hrr(basis_a,basis_b,unique,nz)

L1 = basis_a.L;
L2 = basis_b.L;

Dim1 = (L1+1)*(L1+2)/2;
Dim2 = (L2+1)*(L2+2)/2;
    
tempP = [0;0;0];
tempp = 0;

S = zeros(Dim1,Dim2);
T = zeros(Dim1,Dim2);
S00 = 0;
RAB = -[basis_b.g(1).x0-basis_a.g(1).x0;basis_b.g(1).y0-basis_a.g(1).y0;basis_b.g(1).z0-basis_a.g(1).z0];

%Loop over uncontracted part
        for nba = 1:basis_a.n
            
            for nbb = 1:basis_b.n
                g1 = basis_a.g(nba);
                g2 = basis_b.g(nbb);
                p = g1.alpha + g2.alpha;
                
                Px = (g1.alpha*g1.x0 + g2.alpha*g2.x0)/p;
                Py = (g1.alpha*g1.y0 + g2.alpha*g2.y0)/p;
                Pz = (g1.alpha*g1.z0 + g2.alpha*g2.z0)/p;
                
                [S00ab,~,Pab,p] = one_electron_ss(g1,g2,RAB,Px,Py,Pz);
                S00 = S00 + S00ab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                tempP = tempP + Pab*S00ab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N; %maybe later I need to add contraction and normalization coefficient 
                tempp = tempp + p*S00ab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                
                RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];

                
                [Sab,Tab] = one_electron_highL(g1,g2,L1,L2,RPA,RAB,nz,unique);
                %[Sab,Tab] = one_electron(basis{a}.g(nba),basis{b}.g(nbb));
                
                S = S + Sab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                T = T + Tab*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                
                %Calculation of Pmean and pmean requires an S-type overlap,
                %even if basis functions have angular momentum
                
                

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
        
        Pcell = tempP/S00;
        pcell = tempp/S00;

%Contracted operations
        

end