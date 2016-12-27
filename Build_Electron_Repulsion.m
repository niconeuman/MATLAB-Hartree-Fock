function gabcd = Build_Electron_Repulsion(basis)
nbasis = size(basis,2);
gabcd=zeros(nbasis,nbasis,nbasis,nbasis);

for a = 1:nbasis %loops over the 1st contracted basis function
    for na=1:basis{a}.n %loops over the number of primitives in the 1st contracted basis function
        aa=basis{a}.g(na).alpha;
        for b = 1:nbasis %loops over the second contracted basis function
            for nb=1:basis{b}.n %loops over number of primitives in the 2nd contracted basis function
                ab = basis{b}.g(nb).alpha;
                p = aa + ab;
                Px = (basis{a}.g(na).x0*aa+basis{b}.g(nb).x0*ab)/p;
                Py = (basis{a}.g(na).y0*aa+basis{b}.g(nb).y0*ab)/p;
                Pz = (basis{a}.g(na).z0*aa+basis{b}.g(nb).z0*ab)/p;
                EabX = gprod_1D(basis{a}.g(na).x0,aa,basis{b}.g(nb).x0,ab);
                EabY = gprod_1D(basis{a}.g(na).y0,aa,basis{b}.g(nb).y0,ab);
                EabZ = gprod_1D(basis{a}.g(na).z0,aa,basis{b}.g(nb).z0,ab);
                A_AB = EabX*EabY*EabZ*basis{a}.c(na)*basis{b}.c(nb)*basis{a}.g(na).N*basis{b}.g(nb).N;
                for c = 1:nbasis %loops over 3rd contracted basis function
                    for nc = 1:basis{c}.n; %loop over number of primitives in 3rd contracted basis function 
                    ac=basis{c}.g(nc).alpha;
                    for d=1:nbasis %loops over 4th contracted basis function
                        for nd = 1:basis{d}.n; %loops over number of primitives in 4th contracted basis functions.
                        ad = basis{d}.g(nd).alpha;
                        pp = ac+ad;
                        PPx = (basis{c}.g(nc).x0*aa+basis{d}.g(nd).x0*ab)/pp;
                        PPy = (basis{c}.g(nc).y0*aa+basis{d}.g(nd).y0*ab)/pp;
                        PPz = (basis{c}.g(nc).z0*aa+basis{d}.g(nd).z0*ab)/pp;
                        EcdX = gprod_1D(basis{c}.g(nc).x0,ac,basis{d}.g(nd).x0,ad);
                        EcdY = gprod_1D(basis{c}.g(nc).y0,ac,basis{d}.g(nd).y0,ad);
                        EcdZ = gprod_1D(basis{c}.g(nc).z0,ac,basis{d}.g(nd).z0,ad);
                        A_CD = EcdX*EcdY*EcdZ*basis{c}.c(nc)*basis{d}.c(nd)*basis{c}.g(nc).N*basis{d}.g(nd).N;
                        RPPP2 = (Px-PPx)^2+(Py-PPy)^2+(Pz-PPz)^2;
                        alpha = pp*p/(pp+p);
                        gabcd(a,b,c,d) = gabcd(a,b,c,d)+A_AB*A_CD*Boys(0,alpha*RPPP2)*2*pi^2.5/(p*pp*sqrt(p+pp)); 
                        end
                    end
                    end
                end
            end
        end
    end
end



end