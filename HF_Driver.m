close all
tic
%Boys_Table 400,4000,order 0...10

Boys_Table = loadBoysTable;

nz = cell(1,12);
nz{2} = [1;0;0;1;0;1];
nz{3} = [2;1;1;0;0;0;2;1;0;2];
nz{4} = [3;2;2;1;1;1;0;0;0;0;3;2;1;0;3];
nz{5} = [4;3;3;2;2;2;1;1;1;1;0;0;0;0;0;4;3;2;1;0;4];
nz{6} = [5;4;4;3;3;3;2;2;2;2;1;1;1;1;1;0;0;0;0;0;0;5;4;3;2;1;0;5];
nz{7} = [6;5;5;4;4;4;3;3;3;3;2;2;2;2;2;1;1;1;1;1;1;0;0;0;0;0;0;0;6;5;4;3;2;1;0;6];
nz{8} = [7;6;6;5;5;5;4;4;4;4;3;3;3;3;3;2;2;2;2;2;2;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;7;6;5;4;3;2;1;0;7];
nz{9} = [8;7;7;6;6;6;5;5;5;5;4;4;4;4;4;3;3;3;3;3;3;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;8;7;6;5;4;3;2;1;0;8];
nz{10} = [9;8;8;7;7;7;6;6;6;6;5;5;5;5;5;4;4;4;4;4;4;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;9;8;7;6;5;4;3;2;1;0;9];
nz{11} = [10;9;9;8;8;8;7;7;7;7;6;6;6;6;6;5;5;5;5;5;5;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;10;9;8;7;6;5;4;3;2;1;0;10];
nz{12} = [11;10;10;9;9;9;8;8;8;8;7;7;7;7;7;6;6;6;6;6;6;5;5;5;5;5;5;5;4;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;11;10;9;8;7;6;5;4;3;2;1;0;11];

unique = [2, 4,5, 7,8,9, 11,12,13,14, 16,17,18,19,20, 22,23,24,25,26,27,  29,30,31,32,33,34,35, 37,38,39,40,41,42,43,44, 46,47,48,49,50,51,52,53,54, 56,57,58,59,60,61,62,63,64,65, 67,68,69,70,71,72,73,74,75,76,77, 79,80,81,82,83,84,85,86,87,88,89,90]; %up to [ms|ss]



%   Z = [2 6 1 1 1 1]; %H, He
%  Z = [Z Z];
% %Z = [Z Z Z];
%
%acetilene
% Z = [6 6 1 1];
%% AL = [0.05 0 -0.5
%%       -0.05 0 0.5
%%       0.05 0 -1.4
%%       -0.05 0 1.4]*1.8897;
%Charge = 0;
%AL = [0.094486    0.000000   -0.944863
%-0.094486    0.000000    0.944863
% 0.094486    0.000000   -2.645617
%-0.094486    0.000000    2.645617];



%H2O
% Z = [8 1 1];
% AL = [0.00 0.000 0.116
%       0.000 0.751 -0.465
%       0.000 -0.751 -0.465]*1.8897;
%CO
% Z = [6 8];
% r = 1.10*1.8897;
% theta = 30*pi/180;
% phi = 45*pi/180;
% AL = [0.00 0.00 0
%       sin(theta)*cos(phi)*r sin(theta)*sin(phi)*r cos(theta)*r];

%Methane
 Z = [6 1 1 1 1];
 Charge = 0;
 AL = [-0.000000000000   0.000000000000   0.000000000000
      1.183772681898  -1.183771681898  -1.183771681898
      1.183771681898   1.183772681898   1.183771681898
      -1.183771681898   1.183771681898  -1.183771681898
      -1.183771681898  -1.183771681898   1.183771681898];
%Fictituous Molecule with
%      Z = [6 1 1 1 1];
%      Charge = 0;
      AL = [-0.000000000000   0.000000000000   0.000000000000
           2.183772681898  -2.183771681898  -1.143771681898
           1.183771681898   3.183772681898   1.169771681898
           -1.444771681898   1.883771681898  -1.03771681898
           -1.123771681898  -1.783771681898   1.483771681898];
% %
% % %1 Bohr = 0.5291772108 A;
% %
% AL = [-0.000000000000   0.000000000000   0.000000000000
%     1.183772  -1.183772  -1.183772
%     1.183772   1.183772   1.183772
%     -1.183772   1.183772  -1.183772
%     -1.183772  -1.183772   1.1837712];
% AL = AL*1.5;

% Z = [6];
% Charge = 0;
% AL = [-0.000000000000   0.000000000000   0.000000000000];

%Fraction of molecule_004000
% Z = [8 7 6];
% Charge = 0;
% AL = [-0.2849591644	1.7190324354	0.0525888862
% 0.1014649249	0.4895501827	0.5551454214
% -0.8014553823	-0.4157525006	0.401989349]*1.8897;

%Methane + 2 HeH+
% Z = [6 1 1 1 1 2 1 2 1];
% Charge = 2;
% AL = [-0.000000000000   0.000000000000   0.000000000000
%    1.183771681898  -1.183771681898  -1.183771681898
%    1.183771681898   1.183771681898   1.183771681898
%    -1.183771681898   1.183771681898  -1.183771681898
%    -1.183771681898  -1.183771681898   1.183771681898
%    0   0   4
%    0   0.3 4.8
%    0   0   -4
%    0   -0.3 -4.8];

%CHe2H2
% Z = [6 2 2 1 1];
% AL = [-0.000000000000   0.000000000000   0.000000000000
%     1.183771681898  -1.183771681898  -1.183771681898
%     1.183771681898   1.183771681898   1.183771681898
%     -1.183771681898   1.183771681898  -1.183771681898
%     -1.183771681898  -1.183771681898   1.183771681898];
%
%  Z = [Z Z];
%  AL = [AL;AL+4*ones(size(AL))];

% AL =	[-4 -4 -4
%         -0.000000000000   0.000000000000   0.000000000000
%          1.183771681898  -1.183771681898  -1.183771681898
%          1.183771681898   1.183771681898   1.183771681898
%          -1.183771681898   1.183771681898  -1.183771681898
%          -1.183771681898  -1.183771681898   1.183771681898];
% % AL = AL*1.8897;
%  AL = [AL;AL+3*ones(size(AL))];
%AL = [AL;AL+2*ones(size(AL));AL-2*ones(size(AL))];

%  Z = [1 1];
%  AL = [0 0 0;0 0 0.7408478]*1.8897;
%HeH+ hexamer
% Z = [2 1];
% Z = [Z Z Z Z Z Z];
% AL = [0 0 0;0 0 1.4632];
% AL2 = [0.21 0.5 0;0 0 1.3632];
% AL = [AL;AL+2*ones(size(AL));AL-[3 0 0;3 0 0];AL+4*ones(size(AL));AL2;AL2-5*ones(size(AL2))];

Ndistances = 1;
rmsd_gabcd_k = zeros(1,Ndistances);
for ka = 1:Ndistances
    %Z = [2 1];
   %Z = [Z Z Z Z Z Z];
    %AL = [0 0 0;0 0 1.4632];
    %AL = [0 0 0;0 0 k/2]/0.529177;
    %AL2 = [0.21 0.5 0;0 0 1.3632];
   %AL = [AL;AL+2*ones(size(AL));AL-[3 0 0;3 0 0];AL+4*ones(size(AL));AL2;AL2-5*ones(size(AL2))];

[basis,Nel] = Build_Basis(Z,AL,'STO-3G');
Nel = Nel-Charge;

[basis2,~] = Build_Basis_2(Z,AL,'STO-3G');

[Shell_Doublets,NShell_Doublets] = Build_Shell_Doublets(basis);
Shell_List = Build_Shell_List(basis);

pair_data = basis_products(basis,Shell_Doublets);
[pair_data2,RABdata] = basis_products_2(basis2,Shell_Doublets);

[S,T,Pcell,pcell,S00cell] = Build_One_Electron_3(basis,Shell_Doublets,NShell_Doublets);

%Nuclear_Attraction = Build_Nuclear_Attraction_2(basis,Shell_Doublets,NShell_Doublets,AL,Z,Boys_Table);
Nuclear_Attraction = Build_Nuclear_Attraction_3(basis,Shell_List,AL,Z,Boys_Table);

%[Shells,NShells] = Build_Shells(basis);
%gabcd = Build_ERI(basis,Shells,NShells,Boys_Table,pair_data);
%gabcd = Build_ERI_4(basis,Shell_List,Boys_Table,pair_data);
gabcd = Build_ERI_OS(basis,Shell_List,Boys_Table,pair_data);
[gabcd_fast,int_error] = Build_ERI_fast(basis,Shell_List,Boys_Table,pair_data,S,T,Pcell,pcell,S00cell,gabcd);

% t = 0;
% Ncont = Shell_List(end,1);
% gabcd_list = zeros(nb^4,5);
% for i = 1:Ncont
%     for j = 1:Ncont
%         for k = 1:Ncont
%             for l = 1:Ncont
%                 t = t + 1;
%                 gabcd_list(t,:) = [i j k l gabcd(i,j,k,l)];
%             end
%         end
%     end
% end
gabcd_list = reshape(gabcd,[],1);
gabcd_fast_list = reshape(gabcd_fast,[],1);
int_error_list = gabcd_list-gabcd_fast_list;
rmsd_gabcd = sqrt(sum((gabcd_list-gabcd_fast_list).^2))/length(gabcd_list);

%This calculates the rmsd of the integrals whose error lies below a certain
%threshold.
%This is because some of the integrals have to big errors that dominate the
%rmsd calculation.
int_error_list_below_threshold = int_error_list(abs(int_error_list)<0.01);
rmsd_below_threshold = sqrt(sum(int_error_list_below_threshold.^2))/length(int_error_list_below_threshold);
rmsd_gabcd_k(ka) = rmsd_below_threshold;

%
% if rmsd_gabcd < 0.05
%     rmsd_gabcd_k(k) = rmsd_gabcd;
% else
%     rmsd_gabcd_k(k) = rmsd_gabcd_k(k-1)
% end
%gabcd = Build_Electron_Repulsion(basis);

%This calculates the number of integrals which are below a certain
%threshold
threshold = 0.002;
Int_below_threshold = sum(abs(int_error_list)<threshold);
Percent_below_threshold = Int_below_threshold/length(int_error_list)*100;

H0 = T+Nuclear_Attraction;

[Min_Energy,E,ncycle,D,Dinit,G,epsilon,F,Fprime] = SCF(H0,T,Nuclear_Attraction,gabcd,S,Nel,Shell_List);
%Min_Energy = Min_Energy + Nuclear_Repulsion(Z,AL);
E = E + Nuclear_Repulsion(Z,AL);
toc

if ka<4
 figure;
 plot((1:ncycle)',E,'.k','MarkerFaceColor','k');
 xlabel('Cycle');
 ylabel('Energy (Eh)');

Ncont = Shell_List(end,2);
figure;
plot((1:Ncont^4)',gabcd_list,'.k',(1:Ncont^4)',gabcd_fast_list,'-r');% hold on
%plot((1:Ncont^4)',int_error_list,'-b'); %hold on
%plot((1:Ncont^4)',gabcd_fast_list./gabcd_list,'-m');
legend(num2str(rmsd_gabcd_k(ka)));
end

end
% figure
% plot((1:Ndistances)/2,rmsd_gabcd_k,'sqk','MarkerFaceColor','k');
% diff = zeros(912,1);
% for t = 1:912
% diff(t) = (gabcd(gabcd_sherrill(t,1),gabcd_sherrill(t,2),gabcd_sherrill(t,3),gabcd_sherrill(t,4))-gabcd_sherrill(t,5));
% end
% figure; plot(diff);
