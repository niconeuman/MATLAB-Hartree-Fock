function gabcd = contract_ascs(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique)
%June 21st 2016
%This function contracts integrals of the type (LaS|LcS)

% L1 = basis_a.L;
% L2 = basis_b.L;
% L3 = basis_c.L;
% L4 = basis_d.L;
Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
Dim2 = (L2+1)*(L2+2)/2;
Dim3 = (L3+1)*(L3+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
Dim4 = (L4+1)*(L4+2)/2;
Lmax = L1+L2+L3+L4;
order = 0;
gabcd = zeros(Dim1,Dim3,Dim2,Dim4);

    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        temp = zeros(Dim1,Dim3,Dim2,Dim4); %I permute dimensions 2 and 3 just for acsc
               
            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;
                
                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;
                
                RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
%                 RABx = g1.x0-g2.x0;
%                 RABy = g1.y0-g2.y0;
%                 RABz = g1.z0-g2.z0;
                
                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                
                temp2 = zeros(Dim1,Dim3,Dim2,Dim4);
                
                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;
                        
                        temp3 = zeros(Dim1,Dim3,Dim2,Dim4);
                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
                            
                            g4 = basis_d.g(nd);
                            ad = g4.alpha;
                            c4 = basis_d.c(nd);
                            N4 = g4.N;
                            
                            q = ac+ad;
                            Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            Qz = (ac*g3.z0 + ad*g4.z0)/q;
%                             EcdX = gprod_1D(basis_c.g(nc).x0,ac,basis_d.g(nd).x0,ad);
%                             EcdY = gprod_1D(basis_c.g(nc).y0,ac,basis_d.g(nd).y0,ad);
%                             EcdZ = gprod_1D(basis_c.g(nc).z0,ac,basis_d.g(nd).z0,ad);
%                             A_CD = EcdX*EcdY*EcdZ*basis_c.c(nc)*basis_d.c(nd)*basis_c.g(nc).N*basis_d.g(nd).N;
                            RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
%                             RCDx = g3.x0-g4.x0;
%                             RCDy = g3.y0-g4.y0;
%                             RCDz = g3.z0-g4.z0;
                            Wx = (p*Px+q*Qx)/(p+q);
                            Wy = (p*Py+q*Qy)/(p+q);
                            Wz = (p*Pz+q*Qz)/(p+q);
                            
                            RWP = [Wx-Px;Wy-Py;Wz-Pz];
                            RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];
                            
                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
                            
                            
                            
                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);
                            
                            %gabcd(a,b,c,d) = gabcd(a,b,c,d)+A_AB*A_CD*Boys(0,alpha*RPQ2)*2*pi^2.5/(p*pp*sqrt(p+pp)); 
                            
                            temp3 = temp3 + ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz,unique)*c4*N4;
                            
                        end
                        temp2 = temp2 + temp3*c3*N3;
                    end
                    temp = temp + temp2*c2*N2;
            end
            gabcd = gabcd + temp*c1*N1;
        
    end



end