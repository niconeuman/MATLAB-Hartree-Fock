function VAB = electron_nuclear(g1,g2,Cx,Cy,Cz,Z)

a = g1.alpha;
b = g2.alpha;
p = a+b;
P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
C = [Cx,Cy,Cz];
%Johns calls the C vector A, but in Helgaker's book it is C, so in order to
%avoid confusion I will use C (for the coordinates of the nuclei)

RPC2 = sum((C-P).^2);

if (g1.L == 0 && g2.L == 0)

Ex = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
Ey = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
Ez = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);

%This equation has the same form as eq. 9.10.3 of Helgaker's book, provided
%that K_ab^xyz = 1
Kab = Ex*Ey*Ez;
VAB = -1*Z*(2*pi/p)*Interpolated_Boys(0,p*RPC2)*Kab;

%I can now use formulas from Helgaker's book or Obara and Saika's paper
%NO EXCUSES (5 jan 2016)

else
    %recursive_electron_nuclear is a recursive function which outputs
    %matrices of different sizes depending on g1.L and g2.L
    
    a = g1.alpha;
    b = g2.alpha;
    p = a+b;
    P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
    C = [Cx,Cy,Cz];
    RPC2 = sum((C-P).^2);
     
    %Data I need: p, XPA, YPA, ZPA, XPB, YPB, ZPB, RPA2, XPC, YPC, ZPC,Z
    
    %[Ex,p,q,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    %[Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    %[Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %Row vector
    RPA = (P-[g1.x0,g1.y0,g1.z0]);
    RPB = (P-[g2.x0,g2.y0,g2.z0]);
    RPC = P-C;
    %Order of the integrals. it is zero for the final integral, and
    %increases for intermediate integrals
    order = 0;
    
    VAB = recursive_electron_nuclear(a,b,RPA,RPB,RPC,p,RPC2,Z,g1.L,g2.L,order);


end

end