
Z = [1 1]; %H, He
AL = [0 0 0;0 0 1.4];

[basis,N] = Build_Basis(Z,AL);
[S,T] = Build_Overlap(basis);
Nuclear_Attraction = Build_Nuclear_Attraction(basis,AL,Z);
gabcd = Build_Electron_Repulsion(basis);

H0 = T+Nuclear_Attraction;

[Min_Energy,E,ncycle,D] = SCF(H0,gabcd,S,N);
%Min_Energy = Min_Energy + Nuclear_Repulsion(Z,AL);
E = E + Nuclear_Repulsion(Z,AL);
