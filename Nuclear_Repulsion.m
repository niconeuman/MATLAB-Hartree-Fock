function out = Nuclear_Repulsion(Z,AL)
nAtoms = size(Z,2);
Nuc_Rep = 0;
for n = 1:nAtoms
    for m = (n+1):nAtoms
        Nuc_Rep = Nuc_Rep + Z(n)*Z(m)/norm(AL(n,:)-AL(m,:));
    end
end
out = Nuc_Rep;
end