function [S,T,Pcell,pcell] = Build_One_Electron_3(basis,Shell_Doublets,NShell_Doublets)

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
nb = size(basis,1);
Pcell = cell(nb,nb);
pcell = cell(nb,nb);
S00cell = cell(nb,nb);

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
        [Sout,Tout,Pout,pout] = contract_one_electron(basis{a},basis{b},nz,unique);
        
        S(mu_begin:mu_end,nu_begin:nu_end) = Sout;
        T(mu_begin:mu_end,nu_begin:nu_end) = Tout;
        Pcell{a,b} = Pout;
        pcell{a,b} = pout;
        
    elseif (basis{a}.L > 0 && basis{b}.L == 0)
        
        [Sout,Tout,Pout,pout,S00] = contract_one_electron_highL(basis{a},basis{b},unique,nz);
        S(mu_begin:mu_end,nu_begin:nu_end) = Sout;
        T(mu_begin:mu_end,nu_begin:nu_end) = Tout;
        Pcell{a,b} = Pout;
        pcell{a,b} = pout;
        S00cell{a,b} = S00;
        
    elseif (basis{a}.L == 0 && basis{b}.L > 0)
        
        [Sout,Tout,Pout,pout,S00] = contract_one_electron_highL(basis{b},basis{a},unique,nz);
        S(mu_begin:mu_end,nu_begin:nu_end) = Sout';
        T(mu_begin:mu_end,nu_begin:nu_end) = Tout';
        Pcell{a,b} = Pout;
        pcell{a,b} = pout;
        S00cell{a,b} = S00;
        
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
        %Alternative to check the hrr (gives the same result) 10 nov 2017
%         RAB = [basis{b}.g(1).x0-basis{a}.g(1).x0;basis{b}.g(1).y0-basis{a}.g(1).y0;basis{b}.g(1).z0-basis{a}.g(1).z0];
%         
%         [Sabm1,~] = contract_one_electron_highL_hrr(basis{a},basis{b},L1,L2-1,unique,nz);
%         [Sap1bm1,~] = contract_one_electron_highL_hrr(basis{a},basis{b},L1+1,L2-1,unique,nz);
%    
%         ExpSabm1 = int_expand(Sabm1,RAB,L2-1,2);
%             %termS1 = int_reshape(ExpSabm1,L1+1,L2-1,1,unique);
%         termS2 = int_reshape(Sap1bm1,L1+1,L2-1,1,unique);
%         
%         S(mu_begin:mu_end,nu_begin:nu_end) = ExpSabm1+termS2;
        
        if (basis{b}.L > basis{a}.L)
            [Sout,Tout,Pout,pout,S00] = contract_one_electron_hrr(basis{b},basis{a},unique,nz);
            S(mu_begin:mu_end,nu_begin:nu_end) = Sout';
            T(mu_begin:mu_end,nu_begin:nu_end) = Tout';
            
        else
            [Sout,Tout,Pout,pout,S00] = contract_one_electron_hrr(basis{a},basis{b},unique,nz);
            S(mu_begin:mu_end,nu_begin:nu_end) = Sout;
            T(mu_begin:mu_end,nu_begin:nu_end) = Tout;
        end
            Pcell{a,b} = Pout;
            pcell{a,b} = pout;
            S00cell{a,b} = S00;
    end
    
    

end

end