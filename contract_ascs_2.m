function [gascs,gap1scm1s,gam1scm1s,gascm1s,gascm2s] = contract_ascs_2(basis_a,basis_b,basis_c,basis_d,L1,L2,L3,L4,Boys_Table,nz,unique)
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
%Lmax = L1+L2+L3+L4;
order = 0;
gascs = zeros(Dim1,Dim3);
gap1scm1s = zeros((L1+2)*(L1+3)/2,(L3)*(L3+1)/2);
gam1scm1s = zeros((L1)*(L1+1)/2,(L3)*(L3+1)/2);
gascm1s = zeros((L1+1)*(L1+2)/2,(L3)*(L3+1)/2);
if L3 > 1
gascm2s = zeros((L1+1)*(L1+2)/2,(L3-1)*(L3)/2);
else
gascm2s = [];
end

    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
        %temp = zeros(Dim1,Dim3,Dim2,Dim4); %I permute dimensions 2 and 3 just for ascs
        gascs_temp = zeros(Dim1,Dim3);
        gap1scm1s_temp = zeros((L1+2)*(L1+3)/2,(L3)*(L3+1)/2);
        gam1scm1s_temp = zeros((L1)*(L1+1)/2,(L3)*(L3+1)/2);
        gascm1s_temp = zeros((L1+1)*(L1+2)/2,(L3)*(L3+1)/2);
        if L3 > 1
        gascm2s_temp = zeros((L1+1)*(L1+2)/2,(L3-1)*(L3)/2);
        else
        gascm2s_temp = [];
        end
        
        
            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                g2 = basis_b.g(nb);
                ab = g2.alpha;
                c2 = basis_b.c(nb);
                N2 = g2.N;
                
                p = aa + ab;
                Px = (aa*g1.x0 + ab*g2.x0)/p;
                Py = (aa*g1.y0 + ab*g2.y0)/p;
                Pz = (aa*g1.z0 + ab*g2.z0)/p;
                
                %RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
                RPA = -ab/p*RAB;
%                 RABx = g1.x0-g2.x0;
%                 RABy = g1.y0-g2.y0;
%                 RABz = g1.z0-g2.z0;
                
                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                
                %temp2 = zeros(Dim1,Dim3,Dim2,Dim4);
                gascs_temp2 = zeros(Dim1,Dim3);
                gap1scm1s_temp2 = zeros((L1+2)*(L1+3)/2,(L3)*(L3+1)/2);
                gam1scm1s_temp2 = zeros((L1)*(L1+1)/2,(L3)*(L3+1)/2);
                gascm1s_temp2 = zeros((L1+1)*(L1+2)/2,(L3)*(L3+1)/2);
                if L3 > 1
                gascm2s_temp2 = zeros((L1+1)*(L1+2)/2,(L3-1)*(L3)/2);
                else
                gascm2s_temp2 = [];
                end                
                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;
                        
                        %temp3 = zeros(Dim1,Dim3,Dim2,Dim4);
                        gascs_temp3 = zeros(Dim1,Dim3);
                        gap1scm1s_temp3 = zeros((L1+2)*(L1+3)/2,(L3)*(L3+1)/2);
                        gam1scm1s_temp3 = zeros((L1)*(L1+1)/2,(L3)*(L3+1)/2);
                        gascm1s_temp3 = zeros((L1+1)*(L1+2)/2,(L3)*(L3+1)/2);
                        if L3 > 1
                        gascm2s_temp3 = zeros((L1+1)*(L1+2)/2,(L3-1)*(L3)/2);
                        else
                        gascm2s_temp3 = [];
                        end        
                        
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
                            %RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                            RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
                            RQC = -ad/q*RCD;
                            %Qx-Cx = (c*Cx+d*Dx)/q-q*Cx/q =
                            %c*Cx/q+d*Dx/q-q*Cx/q =
                            %(c-c-d)*Cx/q+d*Dx/q=(-d)*Cx/q+d*Dx/q =
                            %-d/q*(Cx-Dx)
                            
%                             RCDx = g3.x0-g4.x0;
%                             RCDy = g3.y0-g4.y0;
%                             RCDz = g3.z0-g4.z0;
%                             Wx = (p*Px+q*Qx)/(p+q);
%                             Wy = (p*Py+q*Qy)/(p+q);
%                             Wz = (p*Pz+q*Qz)/(p+q);
                            
                            %RWP = [Wx-Px;Wy-Py;Wz-Pz];
                            %RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];
                            
                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
                            
                            
                            
                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RWP = -q/(p+q)*RPQ;
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);
                            
% Kab = vpa(Kab,20);
% Kcd = vpa(Kcd,20);
% RPA = vpa(RPA,20);
% RQC = vpa(RQC,20);
% RWP = vpa(RWP,20);
% p = vpa(p,20);
% q = vpa(q,20);
% alpha = vpa(alpha,20);
% RPQ2 = vpa(RPQ2,20);

                            
                            %temp3 = temp3 + ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz,unique)*c4*N4;
                            %temp3 = temp3 + ascs_3minimal(Kab,Kcd,RPA,RQC,RWP,p,q,alpha,RPQ2,L1,L3,order,Boys_Table,nz,unique)*c4*N4;
                            [gascs_temp4,gap1scm1s_temp4,gam1scm1s_temp4,gascm1s_temp4,gascm2s_temp4] = ascs_4(Kab,Kcd,RPA,RQC,RWP,p,q,alpha,RPQ2,L1,L3,order,Boys_Table,nz,unique);
                            
                            gascs_temp3 = gascs_temp3+gascs_temp4*c4*N4;
                            gap1scm1s_temp3 = gap1scm1s_temp3+gap1scm1s_temp4*c4*N4;
                            gam1scm1s_temp3 = gam1scm1s_temp3+gam1scm1s_temp4*c4*N4;
                            gascm1s_temp3 = gascm1s_temp3+gascm1s_temp4*c4*N4;
                            %Now I don't need the if. In any case I'm
                            %multiplying an empty matrix
                            gascm2s_temp3 = gascm2s_temp3+gascm2s_temp4*c4*N4;

                        
                        end
                        %temp2 = temp2 + temp3*c3*N3;
                            gascs_temp2 = gascs_temp2+gascs_temp3*c3*N3;
                            gap1scm1s_temp2 = gap1scm1s_temp2+gap1scm1s_temp3*c3*N3;
                            gam1scm1s_temp2 = gam1scm1s_temp2+gam1scm1s_temp3*c3*N3;
                            gascm1s_temp2 = gascm1s_temp2+gascm1s_temp3*c3*N3;
                            gascm2s_temp2 = gascm2s_temp2+gascm2s_temp3*c3*N3;
                    end
                    %temp = temp + temp2*c2*N2;
                        gascs_temp = gascs_temp+gascs_temp2*c2*N2;
                        gap1scm1s_temp = gap1scm1s_temp+gap1scm1s_temp2*c2*N2;
                        gam1scm1s_temp = gam1scm1s_temp+gam1scm1s_temp2*c2*N2;
                        gascm1s_temp = gascm1s_temp+gascm1s_temp2*c2*N2;
                        gascm2s_temp = gascm2s_temp+gascm2s_temp2*c2*N2;
            end
%            gabcd = gabcd + temp*c1*N1;
                    gascs = gascs+gascs_temp*c1*N1;
                    gap1scm1s = gap1scm1s+gap1scm1s_temp*c1*N1;
                    gam1scm1s = gam1scm1s+gam1scm1s_temp*c1*N1;
                    gascm1s = gascm1s+gascm1s_temp*c1*N1;
                    gascm2s = gascm2s+gascm2s_temp*c1*N1;        
    end



end