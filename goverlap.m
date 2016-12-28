function S = goverlap(g1,g2)
%I need to modify this function to include overlap for functions with
%angular momentum larger than zero

if (g1.L == 0 && g2.L == 0)
    Ex = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
    Ey = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    Ez = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    S = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3*g1.N*g2.N;

elseif (g1.L == 1 && g2.L ==0)
    [Ex,p,q,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    [Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %I'm explicitly writing the recursion relations, because there are not
    %so many cases.
    S00 = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3*g1.N*g2.N;
    Sx0 = (Px-g1.x0)*S00; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    Sy0 = (Py-g1.y0)*S00;
    Sz0 = (Pz-g1.z0)*S00;
    S = [Sx0;Sy0;Sz0];
elseif (g1.L == 0 && g2.L == 1)
    [Ex,p,q,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    [Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %I'm explicitly writing the recursion relations, because there are not
    %so many cases.
    S00 = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3*g1.N*g2.N;
    S0x = (Px-g2.x0)*S00; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    S0y = (Py-g2.y0)*S00;
    S0z = (Pz-g2.z0)*S00;
    S = [S0x,S0y,S0z]; %row vector
elseif (g1.L == 1 && g2.L == 1)
    [Ex,p,q,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
    [Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    [Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    %I'm explicitly writing the recursion relations, because there are not
    %so many cases.
    S00 = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3*g1.N*g2.N;
    
    Sx0 = (Px-g1.x0)*S00; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    Sy0 = (Py-g1.y0)*S00;
    Sz0 = (Pz-g1.z0)*S00;
    
    S0x = (Px-g2.x0)*S00; %Px,Py,Pz depend on alphas, so this has to be performed on primitives
    S0y = (Py-g2.y0)*S00;
    S0z = (Pz-g2.z0)*S00;
    
    Sxx = (Px-g1.x0)*S0x+1/2/p*(1*Sx0);
    Syx = (Py-g1.y0)*S0x+1/2/p*(1*Sy0);
    Szx = (Pz-g1.z0)*S0x+1/2/p*(1*Sz0);
    
    Sxy = (Px-g1.x0)*S0y+1/2/p*(1*Sx0);
    Syy = (Py-g1.y0)*S0y+1/2/p*(1*Sy0);
    Szy = (Pz-g1.z0)*S0y+1/2/p*(1*Sz0);
    
    Sxz = (Px-g1.x0)*S0z+1/2/p*(1*Sx0);
    Syz = (Py-g1.y0)*S0z+1/2/p*(1*Sy0);
    Szz = (Pz-g1.z0)*S0z+1/2/p*(1*Sz0);
    
    S = [Sxx,Sxy,Sxz;Syx,Syy,Syz;Szx,Szy,Szz]; %row vector
end
end