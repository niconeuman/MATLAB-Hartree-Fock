function [S,T,P,p] = one_electron_2(g1,g2)
%December 30th 2016
%Is it worth to use recursion relations for the one electron integrals?
%If I have f-functions (L = 3), I need the cases
%0,0 -------------------------------1 case
%0,1 1,1 1,0 -----------------------3 cases (4 total)
%2,0 2,1 2,2 1,2 0,2 ---------------5 cases (9 total)
%3,0 3,1 3,2 3,3 2,3 1,3 0,3 -------7 cases (16 total)

%I can certainly program all these cases, but when I add g-functions there
%will be 25 cases. And each case will contain all previous data.

%My choice is that one_electron will only deal with s and p functions
%and one_electron_highL will deal recursively with higher angular momentum
%functions

a = g1.alpha;
b = g2.alpha;
p = a+b;
%P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
R12sq = (g1.x0-g2.x0)^2+(g1.y0-g2.y0)^2+(g1.z0-g2.z0)^2;

if (g1.L == 0 && g2.L == 0)
    [Ex,~,~,Px,~] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
    [Ey,~,~,Py,~] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,~] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    
    S = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3;
    T = S*(a*b/p*(3-2*a*b/p*R12sq));
    P = [Px;Py;Pz];
elseif (g1.L == 1 && g2.L ==0)
    [Ex,~,~,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    [Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %I'm explicitly writing the recursion relations, because there are not
    %so many cases.
    pp = sqrt(pi/(g1.alpha + g2.alpha));
    
    S00x = Ex*pp;
    S00y = Ey*pp;
    S00z = Ez*pp;
    %S00 = S00x*S00y*S00z;
    T00x = S00x*(a*b/p*(1-2*a*b/p*Qx^2));
    T00y = S00y*(a*b/p*(1-2*a*b/p*Qy^2));
    T00z = S00z*(a*b/p*(1-2*a*b/p*Qz^2));
    
    Sx0 = -(Px-g1.x0)*S00x; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    Sy0 = -(Py-g1.y0)*S00y;
    Sz0 = -(Pz-g1.z0)*S00z;
    S = [Sx0*S00y*S00z;
        S00x*Sy0*S00z;
        S00x*S00y*Sz0];
    
    Tx0 = -(Px-g1.x0)*T00x+b/p*(2*a*Sx0);
    Ty0 = -(Py-g1.y0)*T00y+b/p*(2*a*Sy0);
    Tz0 = -(Pz-g1.z0)*T00z+b/p*(2*a*Sz0);
    
    Tp_x0 = Tx0*S00y*S00z+Sx0*T00y*S00z+Sx0*S00y*T00z;
    Tp_y0 = T00x*Sy0*S00z+S00x*Ty0*S00z+S00x*Sy0*T00z;
    Tp_z0 = T00x*S00y*Sz0+S00x*T00y*Sz0+S00x*S00y*Tz0;
    
    T = [Tp_x0;Tp_y0;Tp_z0];
    P = [Px;Py;Pz];
    
elseif (g1.L == 0 && g2.L == 1)
    [Ex,~,~,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    [Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %I'm explicitly writing the recursion relations, because there are not
    %so many cases.
    pp = sqrt(pi/(g1.alpha + g2.alpha));
    
    S00x = Ex*pp;
    S00y = Ey*pp;
    S00z = Ez*pp;
    %S00 = S00x*S00y*S00z;
    T00x = S00x*(a*b/p*(1-2*a*b/p*Qx^2));
    T00y = S00y*(a*b/p*(1-2*a*b/p*Qy^2));
    T00z = S00z*(a*b/p*(1-2*a*b/p*Qz^2));
    
    S0x = -(Px-g2.x0)*S00x; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    S0y = -(Py-g2.y0)*S00y;
    S0z = -(Pz-g2.z0)*S00z;
    S = [S0x*S00y*S00z,S00x*S0y*S00z,S00x*S00y*S0z]; %row vector
    
    T0x = -(Px-g2.x0)*T00x+a/p*(2*b*S0x);
    T0y = -(Py-g2.y0)*T00y+a/p*(2*b*S0y);
    T0z = -(Pz-g2.z0)*T00z+a/p*(2*b*S0z);
    
    T0p_x = T0x*S00y*S00z+S0x*T00y*S00z+S0x*S00y*T00z;
    T0p_y = T00x*S0y*S00z+S00x*T0y*S00z+S00x*S0y*T00z;
    T0p_z = T00x*S00y*S0z+S00x*T00y*S0z+S00x*S00y*T0z;
    
    T = [T0p_x,T0p_y,T0p_z];
    P = [Px;Py;Pz];
elseif (g1.L == 1 && g2.L == 1)
    [Ex,~,~,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    [Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %I'm explicitly writing the recursion relations, because there are not
    %so many cases.
    pp = sqrt(pi/(g1.alpha + g2.alpha));
    
    S00x = Ex*pp;
    S00y = Ey*pp;
    S00z = Ez*pp;
    %S00 = S00x*S00y*S00z;
    
    T00x = S00x*(a*b/p*(1-2*a*b/p*Qx^2));
    T00y = S00y*(a*b/p*(1-2*a*b/p*Qy^2));
    T00z = S00z*(a*b/p*(1-2*a*b/p*Qz^2));
    
    Sx0 = -(Px-g1.x0)*S00x; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    Sy0 = -(Py-g1.y0)*S00y;
    Sz0 = -(Pz-g1.z0)*S00z;
    
    S0x = -(Px-g2.x0)*S00x; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    S0y = -(Py-g2.y0)*S00y;
    S0z = -(Pz-g2.z0)*S00z;
    
    %Sign of Px-g1.x0, etc should be reversed, according to Helgaker's
    %convention
    Sxx = -(Px-g2.x0)*Sx0+1/2/p*(1*S00x);
    Syy = -(Py-g2.y0)*S0y+1/2/p*(1*S00y);
    Szz = -(Pz-g2.z0)*S0z+1/2/p*(1*S00z);
    
    Spx_px = Sxx*S00y*S00z;
    Spx_py = Sx0*S0y*S00z;
    Spx_pz = Sx0*S00y*S0z;
    
    Spy_px = S0x*Sy0*S00z;
    Spy_py = S00x*Syy*S00z;
    Spy_pz = S00x*Sy0*S0z;
    
    Spz_px = S0x*S00y*Sz0;
    Spz_py = S00x*S0y*Sz0;
    Spz_pz = S00x*S00y*Szz;
    
    S = [Spx_px,Spx_py,Spx_pz;
         Spy_px,Spy_py,Spy_pz;
         Spz_px,Spz_py,Spz_pz];
    
    %This seems to be OK, according to my writing down of Helgaker's relations 
    Tx0 = -(Px-g1.x0)*T00x+b/p*(2*a*Sx0);
    Ty0 = -(Py-g1.y0)*T00y+b/p*(2*a*Sy0);
    Tz0 = -(Pz-g1.z0)*T00z+b/p*(2*a*Sz0);
    
    T0x = -(Px-g2.x0)*T00x+a/p*(2*b*S0x);
    T0y = -(Py-g2.y0)*T00y+a/p*(2*b*S0y);
    T0z = -(Pz-g2.z0)*T00z+a/p*(2*b*S0z);
    
    Txx = -(Px-g2.x0)*Tx0+1/2/p*(1*T00x+0)+a/p*(2*b*Sxx-0);
    Tyy = -(Py-g2.y0)*Ty0+1/2/p*(1*T00y+0)+a/p*(2*b*Syy-0);
    Tzz = -(Pz-g2.z0)*Tz0+1/2/p*(1*T00z+0)+a/p*(2*b*Szz-0);
    %up to here, ok
    
    Tp_xp_x = Txx*S00y*S00z+Sxx*T00y*S00z+Sxx*S00y*T00z; %i,j = 1, others zero
    Tp_xp_y = Tx0*S0y*S00z+Sx0*T0y*S00z+Sx0*S0y*T00z; %i,l = 1, others zero   
    Tp_xp_z = Tx0*S00y*S0z+Sx0*T00y*S0z+Sx0*S00y*T0z; %i,n = 1, others zero
    Tp_yp_x = T0x*Sy0*S00z+S0x*Ty0*S00z+S0x*Sy0*T00z; %k,j = 1, others zero
    Tp_yp_y = T00x*Syy*S00z+S00x*Tyy*S00z+S00x*Syy*T00z; %k,l = 1, others zero
    Tp_yp_z = T00x*Sy0*S0z+S00x*Ty0*S0z+S00x*Sy0*T0z; %k, n = 1
    
    Tp_zp_x = T0x*S00y*Sz0+S0x*T00y*Sz0+S0x*S00y*Tz0;%m,j = 1
    Tp_zp_y = T00x*S0y*Sz0+S00x*T0y*Sz0+S00x*S0y*Tz0;%m, l = 1
    Tp_zp_z = T00x*S00y*Szz+S00x*T00y*Szz+S00x*S00y*Tzz;%m,n = 1
    
    
    
    T = [Tp_xp_x,Tp_xp_y,Tp_xp_z;Tp_yp_x,Tp_yp_y,Tp_yp_z;Tp_zp_x,Tp_zp_y,Tp_zp_z];
    P = [Px;Py;Pz];

end