function VAB = electron_nuclear_2(g1,g2,Cx,Cy,Cz,Z,Boys_Table)

a = g1.alpha;
b = g2.alpha;
p = a+b;
%P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
Px = (a*g1.x0 + b*g2.x0)/p;
Py = (a*g1.y0 + b*g2.y0)/p;
Pz = (a*g1.z0 + b*g2.z0)/p;
%C = [Cx,Cy,Cz];
%Johns calls the C vector A, but in Helgaker's book it is C, so in order to
%avoid confusion I will use C (for the coordinates of the nuclei)

RPC2 = (Px-Cx)^2+(Py-Cy)^2+(Pz-Cz)^2;


%Ex = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
%Ey = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
%Ez = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);

Kab = gprod(g1.x0,g1.y0,g1.z0,g1.alpha,g2.x0,g2.y0,g2.z0,g2.alpha);

%This equation has the same form as eq. 9.10.3 of Helgaker's book, provided
%that K_ab^xyz = 1
%Kab = Ex*Ey*Ez;
VAB = -1*Z*(2*pi/p)*Interpolated_Boys(0,p*RPC2,Boys_Table)*Kab;
%VAB = -1*Z*(2*pi/p)*Boys(0,p*RPC2)*Kab;


end