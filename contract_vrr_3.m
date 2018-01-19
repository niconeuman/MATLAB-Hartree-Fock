function [gLSSS,gLm1SSS,xValues,DxValues,indexValues,boysValues] = contract_vrr_3(basis_a,basis_b,basis_c,basis_d,L1,Boys_Table) %nz was replaced by ~

Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc

gLSSS = zeros((L1+1)*(L1+2)/2,1);
gLm1SSS = zeros((L1)*(L1+1)/2,1);
%fun is a function handle. I could use VRR, HRR, etc



%The first part of this function loops over all primitives and calculates
%the x_index and Dx values needed for Boys function calculation.
%If I have angular momentum La, I need Boys function up to La+1.
%So
Nints = basis_a.n*basis_b.n*basis_c.n*basis_d.n;
boysValues = zeros(Nints,L1+1);
xValues = zeros(Nints,1);
indexValues = zeros(Nints,1);
xIndexValues = zeros(Nints,1);
DxValues = zeros(Nints,1);
% Dx2Values = zeros(Nints,1);
% Dx3Values = zeros(Nints,1);
% Dx4Values = zeros(Nints,1);
% Dx5Values = zeros(Nints,1);
xstep = 0.1; %Nlast/Npoints
    


t = 1;
    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        g1 = basis_a.g(na);
        aa = g1.alpha;
        c1 = basis_a.c(na);
        N1 = g1.N;
   
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
                %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
%                 RABx = g1.x0-g2.x0;
%                 RABy = g1.y0-g2.y0;
%                 RABz = g1.z0-g2.z0;
                
                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                

                
                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
                        g3 = basis_c.g(nc);
                        ac = g3.alpha;
                        c3 = basis_c.c(nc);
                        N3 = g3.N;
                        

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
%                             RCDx = g3.x0-g4.x0;
%                             RCDy = g3.y0-g4.y0;
%                             RCDz = g3.z0-g4.z0;
                            Wx = (p*Px+q*Qx)/(p+q);
                            Wy = (p*Py+q*Qy)/(p+q);
                            Wz = (p*Pz+q*Qz)/(p+q);
                            
                            RWP = [Wx-Px;Wy-Py;Wz-Pz];
                            %RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];
                            
                            rhoCD = ac*ad/q;
                            Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
                            
                            
                            
                            RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);
                            
                            x = alpha*RPQ2;

                            
                            xValues(t) = x;
                            t = t + 1;
                            %(na-1)*basis_b.n*basis_c.n*basis_d.n+(nb-1)*basis_c.n*basis_d.n+(nc-1)*basis_d.n+nd
                            

                            
                            
                            
                        end

                    end

            end

        
    end


indexValues = floor(xValues/xstep)+1;
xIndexValues = (indexValues-1)*xstep;
DxValues = (xValues-xIndexValues); %Difference which enters the Taylor expansion                            
    
    

Dx2Values = DxValues.*DxValues;
Dx3Values = Dx2Values.*DxValues;
Dx4Values = Dx2Values.*Dx2Values;
Dx5Values = Dx3Values.*Dx2Values;
%  Dx6 = Dx3*Dx3;
%  Dx7 = Dx4*Dx3;
 Dx2_o_2_Values = .5*Dx2Values;
 Dx3_o_6_Values = .166666666666667*Dx3Values;
 Dx4_o_24_Values = 4.166666666666666e-2*Dx4Values;
 Dx5_o_120_Values = 8.333333333333334e-3*Dx5Values;
%  seven_hundred_twentieth_Dx6 = 1.388888888888889e-3*Dx6;
%  Dx7_over_5040 = 1.984126984126984e-04*Dx7;

boysValues(:,1) = Boys_Table(indexValues,1)-Boys_Table(indexValues,2).*DxValues+...
    Boys_Table(indexValues,3).*Dx2_o_2_Values-Boys_Table(indexValues,4).*Dx3_o_6_Values+...
    Boys_Table(indexValues,5).*Dx4_o_24_Values-Boys_Table(indexValues,6).*Dx5_o_120_Values;
boysValues(:,2) = Boys_Table(indexValues,2)-Boys_Table(indexValues,3).*DxValues+...
    Boys_Table(indexValues,4).*Dx2_o_2_Values-Boys_Table(indexValues,5).*Dx3_o_6_Values+...
    Boys_Table(indexValues,6).*Dx4_o_24_Values-Boys_Table(indexValues,7).*Dx5_o_120_Values;




% for k = 1:Nints                               
% boysValues(k,1) = Boys_Table(indexValues(k),1)-Boys_Table(indexValues(k),2)*DxValues(k)+...
%     Boys_Table(indexValues(k),3)*Dx2_o_2_Values(k)-Boys_Table(indexValues(k),4)*Dx3_o_6_Values(k)+...
%     Boys_Table(indexValues(k),5)*Dx4_o_24_Values(k)-Boys_Table(indexValues(k),6)*Dx5_o_120_Values(k);
% boysValues(k,2) = Boys_Table(indexValues(k),2)-Boys_Table(indexValues(k),3)*DxValues(k)+...
%     Boys_Table(indexValues(k),4)*Dx2_o_2_Values(k)-Boys_Table(indexValues(k),5)*Dx3_o_6_Values(k)+...
%     Boys_Table(indexValues(k),6)*Dx4_o_24_Values(k)-Boys_Table(indexValues(k),7)*Dx5_o_120_Values(k);
% 
% end    
    

%     for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
%         g1 = basis_a.g(na);
%         aa = g1.alpha;
%         c1 = basis_a.c(na);
%         N1 = g1.N;
%         gLSSS_temp = zeros((L1+1)*(L1+2)/2,1);
%         gLm1SSS_temp = zeros((L1)*(L1+1)/2,1);       
%             for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
%                 g2 = basis_b.g(nb);
%                 ab = g2.alpha;
%                 c2 = basis_b.c(nb);
%                 N2 = g2.N;
%                 
%                 p = aa + ab;
%                 Px = (aa*g1.x0 + ab*g2.x0)/p;
%                 Py = (aa*g1.y0 + ab*g2.y0)/p;
%                 Pz = (aa*g1.z0 + ab*g2.z0)/p;
%                 
%                 RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
%                 %RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
%                 RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
% %                 RABx = g1.x0-g2.x0;
% %                 RABy = g1.y0-g2.y0;
% %                 RABz = g1.z0-g2.z0;
%                 
%                 rhoAB = aa*ab/p;
%                 Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
%                 
%                 gLSSS_temp2 = zeros((L1+1)*(L1+2)/2,1);
%                 gLm1SSS_temp2 = zeros((L1)*(L1+1)/2,1); 
%                 
%                     for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
%                         g3 = basis_c.g(nc);
%                         ac = g3.alpha;
%                         c3 = basis_c.c(nc);
%                         N3 = g3.N;
%                         
%                         gLSSS_temp3 = zeros((L1+1)*(L1+2)/2,1);
%                         gLm1SSS_temp3 = zeros((L1)*(L1+1)/2,1); 
%                         for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
%                             
%                             g4 = basis_d.g(nd);
%                             ad = g4.alpha;
%                             c4 = basis_d.c(nd);
%                             N4 = g4.N;
%                             
%                             q = ac+ad;
%                             Qx = (ac*g3.x0 + ad*g4.x0)/q;
%                             Qy = (ac*g3.y0 + ad*g4.y0)/q;
%                             Qz = (ac*g3.z0 + ad*g4.z0)/q;
% %                             EcdX = gprod_1D(basis_c.g(nc).x0,ac,basis_d.g(nd).x0,ad);
% %                             EcdY = gprod_1D(basis_c.g(nc).y0,ac,basis_d.g(nd).y0,ad);
% %                             EcdZ = gprod_1D(basis_c.g(nc).z0,ac,basis_d.g(nd).z0,ad);
% %                             A_CD = EcdX*EcdY*EcdZ*basis_c.c(nc)*basis_d.c(nd)*basis_c.g(nc).N*basis_d.g(nd).N;
%                             %RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
%                             %RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
%                             RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
% %                             RCDx = g3.x0-g4.x0;
% %                             RCDy = g3.y0-g4.y0;
% %                             RCDz = g3.z0-g4.z0;
%                             Wx = (p*Px+q*Qx)/(p+q);
%                             Wy = (p*Py+q*Qy)/(p+q);
%                             Wz = (p*Pz+q*Qz)/(p+q);
%                             
%                             RWP = [Wx-Px;Wy-Py;Wz-Pz];
%                             %RWQ = [Wx-Qx;Wy-Qy;Wz-Qz];
%                             
%                             rhoCD = ac*ad/q;
%                             Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
%                             
%                             
%                             
%                             RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
%                             RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
%                             alpha = q*p/(q+p);
%                             
%                             %gabcd(a,b,c,d) = gabcd(a,b,c,d)+A_AB*A_CD*Boys(0,alpha*RPQ2)*2*pi^2.5/(p*pp*sqrt(p+pp)); 
%                             %temp3 = temp3 + vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table)*c4*N4;
%                             %For checking the different recursion
%                             %strategies
%                             %temp3 = temp3 + vrr_1old(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,[0;0;0],[0;0;0],[0;0;0],RWP,[0;0;0],p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,L1,order,Boys_Table,nz)*c4*N4;
%                             %temp3 = temp3 + vrr_1(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)*c4*N4;
%                             
%                             [gLSSS_temp4,~,gLm1SSS_temp4,~,~] = vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table);
%                             gLSSS_temp3 = gLSSS_temp3 + gLSSS_temp4*c4*N4;
%                             gLm1SSS_temp3 = gLm1SSS_temp3 + gLm1SSS_temp4*c4*N4;
%                         end
%                         gLSSS_temp2 = gLSSS_temp2 + gLSSS_temp3*c3*N3;
%                         gLm1SSS_temp2 = gLm1SSS_temp2 + gLm1SSS_temp3*c3*N3;
%                     end
%                     gLSSS_temp = gLSSS_temp + gLSSS_temp2*c2*N2;
%                     gLm1SSS_temp = gLm1SSS_temp + gLm1SSS_temp2*c2*N2;
%             end
%             gLSSS = gLSSS + gLSSS_temp*c1*N1;
%             gLm1SSS = gLm1SSS + gLm1SSS_temp*c1*N1;
%         
%     end







end