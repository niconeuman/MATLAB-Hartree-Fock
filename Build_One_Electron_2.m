function [S,T,Pcell,pcell] = Build_One_Electron_2(basis,Shell_Doublets,NShell_Doublets)

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
%The cell of vectors P will contain the gravity centers of two CGTOs.
%This center will be the sum of P(nba,nbb)*S(nba,nbb)/Sab
Pcell = cell(Ncont,Ncont);
pcell = cell(Ncont,Ncont);
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

L1 = basis{a}.L;
L2 = basis{b}.L;
    if (basis{a}.L == 0 && basis{b}.L == 0)
        tempP = [0;0;0];
        tempp = 0;
        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                [Sab,Tab,Pab,p] = one_electron_2(basis{a}.g(nba),basis{b}.g(nbb));
                %[Sab,Tab] = one_electron(basis{a}.g(nba),basis{b}.g(nbb));
                
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                T(mu_begin:mu_end,nu_begin:nu_end) = T(mu_begin:mu_end,nu_begin:nu_end) + Tab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                tempP = tempP + Pab*Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N; %maybe later I need to add contraction and normalization coefficient 
                tempp = tempp + p*Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
            end
        end
        %For matrix elements larger than (1x1) I need to modify this in
        %some way.
        Pcell{mu_begin:mu_end,nu_begin:nu_end} = tempP/S(mu_begin:mu_end,nu_begin:nu_end);
        pcell{mu_begin:mu_end,nu_begin:nu_end} = tempp/S(mu_begin:mu_end,nu_begin:nu_end);
    elseif (basis{a}.L > 0 && basis{b}.L == 0)
        RAB = [basis{b}.g(1).x0-basis{a}.g(1).x0;basis{b}.g(1).y0-basis{a}.g(1).y0;basis{b}.g(1).z0-basis{a}.g(1).z0];
        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                g1 = basis{a}.g(nba);
                g2 = basis{b}.g(nbb);
                p = g1.alpha + g2.alpha;
                
                Px = (g1.alpha*g1.x0 + g2.alpha*g2.x0)/p;
                Py = (g1.alpha*g1.y0 + g2.alpha*g2.y0)/p;
                Pz = (g1.alpha*g1.z0 + g2.alpha*g2.z0)/p;
                
                RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                [Sab,Tab] = one_electron_highL(g1,g2,L1,L2,RPA,RAB,nz,unique);
                
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                T(mu_begin:mu_end,nu_begin:nu_end) = T(mu_begin:mu_end,nu_begin:nu_end) + Tab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;              
            end
        end
    elseif (basis{a}.L == 0 && basis{b}.L > 0)
        RAB = [basis{b}.g(1).x0-basis{a}.g(1).x0;basis{b}.g(1).y0-basis{a}.g(1).y0;basis{b}.g(1).z0-basis{a}.g(1).z0];
        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                %Permutation of the a and b basis functions, and
                %permutation of the resultant matrices dimensions
                %(transposition)
                %I need to calculate this stuff with basis{b}.g(nbb)
                %permuted with basis{a}.g(nba)
                g1 = basis{a}.g(nba);
                g2 = basis{b}.g(nbb);
                p = g1.alpha + g2.alpha;
                
                Px = (g1.alpha*g1.x0 + g2.alpha*g2.x0)/p;
                Py = (g1.alpha*g1.y0 + g2.alpha*g2.y0)/p;
                Pz = (g1.alpha*g1.z0 + g2.alpha*g2.z0)/p;
                
                %RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                
                [Sab,Tab] = one_electron_highL(g2,g1,L2,L1,RPB,RAB,nz,unique); %The basis functions are permuted, and I use RPB as RPA in the recursion
                Sab = Sab';
                Tab = Tab';
                
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                T(mu_begin:mu_end,nu_begin:nu_end) = T(mu_begin:mu_end,nu_begin:nu_end) + Tab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;              
            end
        end
    else %(basis{a}.L > 0 && basis{b}.L > 0)
        %In here I need to first calculate the contracted (Ls) integrals,
        %and then use a horizontal recursion relation.
        
        %July 4th 2017
        %The horizontal recursion relation for overlap integrals can be
        %simply written for contracted functions: 
        %[a|b+1] = [a|b]*RPA + [a+1|b]
        %However, for kinetic energy integrals, there are two parts:
        %[a|T|b+1] = [a|T|b]*RPA + [a+1|T|b]
        %+2*a*b/p*(([a|b+1]-[a+1|b])-b_i/2/p[a|b-1]+a_i/2/p[a-1|b]
        
        %Sas = zeros((L1+1)*(L1+2)/2,(L2+1)*(L2+2)/2);
        %Tas = zeros((L1+1)*(L1+2)/2,(L2+1)*(L2+2)/2);
        RAB = [basis{b}.g(1).x0-basis{a}.g(1).x0;basis{b}.g(1).y0-basis{a}.g(1).y0;basis{b}.g(1).z0-basis{a}.g(1).z0];
        for nba = 1:basis{a}.n
            for nbb = 1:basis{b}.n
                g1 = basis{a}.g(nba);
                g2 = basis{b}.g(nbb);
                p = g1.alpha + g2.alpha;
                
                Px = (g1.alpha*g1.x0 + g2.alpha*g2.x0)/p;
                Py = (g1.alpha*g1.y0 + g2.alpha*g2.y0)/p;
                Pz = (g1.alpha*g1.z0 + g2.alpha*g2.z0)/p;
                
                RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %RPB = [Px-g2.x0,Py-g2.y0,Pz-g2.z0];
                
                [Sab,Tab] = one_electron_highL(g1,g2,L1,L2,RPA,RAB,nz,unique);
                
                S(mu_begin:mu_end,nu_begin:nu_end) = S(mu_begin:mu_end,nu_begin:nu_end) + Sab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
                T(mu_begin:mu_end,nu_begin:nu_end) = T(mu_begin:mu_end,nu_begin:nu_end) + Tab*basis{a}.c(nba)*basis{b}.c(nbb)*basis{a}.g(nba).N*basis{b}.g(nbb).N;
            end
        end
        %Once I have
        
    end
    
    

end

end