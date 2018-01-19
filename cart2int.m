function InternalCoords = cart2int(CartesianCoords)

%This function reads a cartesian coordinate matrix in the format
% CartesianCoords = [x1 y1 z1;
%                    x2 y2 z2;
%                    x3 y3 z3;
%                    ...      ];

%and generates internal coordinates in the following format
%Atom1 Atom2 Atom3 Atom4 Bond Angle Dihedral
%1
%2     1                 r12
%3     1                 r31
%3     2                 r32
%etc
%3     2     1                <321                %2 being the center atom
%etc

%InternalCoords will be a cell array (for now with only bond distances)

InternalCoords = cell(1);

NAtoms = size(CartesianCoords,1);
%For bond lengths
%I need 3 cols, one for atom1, one for atom2 and one for bondlength
InternalCoords{1} = zeros(NAtoms*(NAtoms-1)/2,3); 
t = 1;
for k = 2:NAtoms
    for l = 1:k-1
        InternalCoords{1}(t,:)=[k  l  norm(CartesianCoords(k,:)-CartesianCoords(l,:))];
        t = t+1;
    end
end

%For bond angles
%I need Atom1 Atom2 Atom3 BondAngle
InternalCoords{2} = zeros((NAtoms*(NAtoms-1)*(NAtoms-2))/6,4); 
t = 1;
for k = 3:NAtoms
    for l = 2:k-1
        for m = 1:l-1
        vkl = CartesianCoords(k,:)-CartesianCoords(l,:);
        vkl = vkl/norm(vkl);
        vlm = CartesianCoords(m,:)-CartesianCoords(l,:);
        vlm = vlm/norm(vlm);
        Angle_klm = acos(vkl*vlm')*180/pi;
        InternalCoords{2}(t,:)=[k  l  m Angle_klm];
        t = t+1;
        end
    end
end





end