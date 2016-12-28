function [basis, N] = Build_Basis(Z,AL)

%Z(N) set of atomic numbers
%AL(N,3) array of atomic coordinates

Natoms = size(Z,2);
N = 0;
nb = 0;

for index = 1:Natoms
    x0 = AL(index,1);
    y0 = AL(index,2);
    z0 = AL(index,3);
    switch Z(index)
        case 1 %Z(index) == 1, H atom, 6-31g, 2 CFG in total
            N = N+1; %number of electrons?
            S = [18.7311370, 0.03349460
                 2.8253937, 0.23472695
                 0.6401217, 0.81375733];
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
            nb = nb + 1;
            basis{nb}.n = 1; %The 6-31g basis set for H has a 3g contracted function and a 1g function
            basis{nb}.c = 1;
            basis{nb}.g = Build_SOrbital(x0,y0,z0,0.1612778);
            basis{nb}.L = 0; %angular momentum type
        case 2 %He atom, 6-31g
                        N = N+2; %number of electrons?
            S = [38.4216340, 0.0237660
                 5.7780300, 0.1546790
                 1.2417740, 0.4696300];
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
            nb = nb + 1;
            basis{nb}.n = 1; %The 6-31g basis set for H has a 3g contracted function and a 1g function
            basis{nb}.c = 1;
            basis{nb}.g = Build_SOrbital(x0,y0,z0,0.2979640);
            basis{nb}.L = 0; %angular momentum type
        case 6 %C atom, 6-31g, 9 CGF in total
            N = N+6; %number of electrons?
            S = [3047.5249000 0.0018347        
                457.3695100 0.0140373        
                103.9486900 0.0688426        
                29.2101550 0.2321844        
                9.2866630 0.4679413        
                3.1639270 0.3623120];
            SP1 = [7.8682724 -0.1193324 0.0689991        
                  1.8812885 -0.1608542 0.3164240        
                  0.5442493 1.1434564 0.7443083]; %3rd column are c coefficients for p
            SP2 = [0.1687144 1.000 1.000];  
            
            nb = nb+1;
            basis{nb}.n = 6; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.g(4) = Build_SOrbital(x0,y0,z0,S(4,1));
            basis{nb}.g(5) = Build_SOrbital(x0,y0,z0,S(5,1));
            basis{nb}.g(6) = Build_SOrbital(x0,y0,z0,S(6,1));
            basis{nb}.L = 0; %angular momentum type
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP1(:,2);
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,SP1(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,SP1(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,SP1(3,1));
            basis{nb}.L = 0; %angular momentum type
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP1(:,3);
            basis{nb}.g(1) = Build_POrbital(x0,y0,z0,SP1(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_POrbital(x0,y0,z0,SP1(2,1));
            basis{nb}.g(3) = Build_POrbital(x0,y0,z0,SP1(3,1));
            basis{nb}.L = 1;
            nb = nb + 1;
            basis{nb}.n = 1;
            basis{nb}.c = SP2(2);
            basis{nb}.g = Build_SOrbital(x0,y0,z0,SP2(1));
            basis{nb}.L = 0;
            nb = nb + 1;
            basis{nb}.n = 1;
            basis{nb}.c = SP2(3);
            basis{nb}.g = Build_POrbital(x0,y0,z0,SP2(1));
            basis{nb}.L = 1;
            
    end
end

end