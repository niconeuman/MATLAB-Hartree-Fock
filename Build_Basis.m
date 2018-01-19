function [basis, N] = Build_Basis(Z,AL,basisset)

%Z(N) set of atomic numbers
%AL(N,3) array of atomic coordinates

Natoms = size(Z,2);
N = 0;
nb = 0;
basis = cell(200,1);
if (strcmp(basisset,'6-31G') || strcmp(basisset,'6-31g'))

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
                 1.2417740, 0.8835300]; %Cambio temporalmente 0.4696300 por otros valores porque S(He1s,He1s) = 0.35
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

elseif strcmp(basisset,'STO-3G')

    for index = 1:Natoms
    x0 = AL(index,1);
    y0 = AL(index,2);
    z0 = AL(index,3);
    switch Z(index)
        case 1 %Z(index) == 1, H atom, 6-31g, 2 CFG in total
            N = N+1; %number of electrons?
            S = [3.42525091,0.15432897       
                 0.62391373,0.53532814       
                 0.16885540,0.44463454];
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
            
        case 2 %He atom, 6-31g
            N = N+2; %number of electrons?
            S = [6.36242139,0.15432897       
                 1.15892300,0.53532814       
                 0.31364979,0.44463454];
            %This is for comparison with Szabo y Ostlund. The difference in
            %the alpha coefficients with those of standard He STO-3G are
            %due to the fact that the authors already correct for the
            %positive charge (mostly in the He atom). This shows that a
            %more flexible basis is needed.
%               S = [3.42525091*1.6875^2,0.15432897       
%                   0.62391373*1.6875^2,0.53532814       
%                   0.16885540*1.6875^2,0.44463454]; 
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
        
        
        case 6 %C atom, 6-31g, 9 CGF in total
            N = N+6; %number of electrons?
            S = [71.6168370, 0.15432897       
                 13.0450960,0.53532814       
                 3.5305122,0.44463454];
            SP = [2.9412494,-0.09996723,0.15591627       
                  0.6834831,0.39951283,0.60768372       
                  0.2222899,0.70011547,0.39195739]; %3rd column are c coefficients for p
            
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
            
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP(:,2);
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,SP(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,SP(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,SP(3,1));
            basis{nb}.L = 0; %angular momentum type
            
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP(:,3);
            basis{nb}.g(1) = Build_POrbital(x0,y0,z0,SP(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_POrbital(x0,y0,z0,SP(2,1));
            basis{nb}.g(3) = Build_POrbital(x0,y0,z0,SP(3,1));
            basis{nb}.L = 1;
            
        case 7 %N atom, STO-3G
            N = N+7; %number of electrons?
            S = [99.1061690,0.15432897       
                 18.0523120,0.53532814       
                 4.8856602,0.44463454];
            SP = [3.7804559,-0.09996723,0.15591627       
                  0.8784966,0.39951283,0.60768372       
                  0.2857144,0.70011547,0.39195739]; %3rd column are c coefficients for p
            
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
            
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP(:,2);
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,SP(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,SP(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,SP(3,1));
            basis{nb}.L = 0; %angular momentum type
            
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP(:,3);
            basis{nb}.g(1) = Build_POrbital(x0,y0,z0,SP(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_POrbital(x0,y0,z0,SP(2,1));
            basis{nb}.g(3) = Build_POrbital(x0,y0,z0,SP(3,1));
            basis{nb}.L = 1;
       case 8 %O atom, STO-3G
            N = N+8; %number of electrons?
            S = [130.7093200,0.15432897       
                 23.8088610,0.53532814       
                 6.4436083,0.44463454];
            SP = [5.0331513,-0.09996723,0.15591627       
                  1.1695961,0.39951283,0.60768372       
                  0.3803890,0.70011547,0.39195739]; %3rd column are c coefficients for p
            
            nb = nb+1;
            basis{nb}.n = 3; %3 basis functions
            basis{nb}.c = S(:,2); %coefficients of each basis function
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,S(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,S(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,S(3,1));
            basis{nb}.L = 0; %angular momentum type
            
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP(:,2);
            basis{nb}.g(1) = Build_SOrbital(x0,y0,z0,SP(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_SOrbital(x0,y0,z0,SP(2,1));
            basis{nb}.g(3) = Build_SOrbital(x0,y0,z0,SP(3,1));
            basis{nb}.L = 0; %angular momentum type
            
            nb = nb + 1;
            basis{nb}.n = 3;
            basis{nb}.c = SP(:,3);
            basis{nb}.g(1) = Build_POrbital(x0,y0,z0,SP(1,1)); %generates all needed data for each basis function in the contracted function
            basis{nb}.g(2) = Build_POrbital(x0,y0,z0,SP(2,1));
            basis{nb}.g(3) = Build_POrbital(x0,y0,z0,SP(3,1));
            basis{nb}.L = 1;
    end
    end
else
    error('Basis set not recognized');
end

basis = basis(1:nb,1);
end