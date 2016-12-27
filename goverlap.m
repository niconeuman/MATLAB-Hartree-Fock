function S = goverlap(g1,g2)
    Ex = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
    Ey = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    Ez = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    S = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3*g1.N*g2.N;


end