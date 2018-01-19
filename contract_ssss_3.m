function gabcd = contract_ssss_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data)

%April 7th 2017
%This function uses pair_data to avoid many computations which are
%redundant. As I am trying to make the code more efficients, for any changes I
%make from contract_ssss to contract_ssss_2, I will leave the old code
%commented above with the %old prefix.

%old
% L1 = basis_a.L;
% L2 = basis_b.L;
% L3 = basis_c.L;
% L4 = basis_d.L;
% Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
% Dim2 = (L2+1)*(L2+2)/2;
% Dim3 = (L3+1)*(L3+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
% Dim4 = (L4+1)*(L4+2)/2;
% Lmax = L1+L2+L3+L4;
% order = 0;
% gabcd = zeros(Dim1,Dim2,Dim3,Dim4);
%end old

%As it is a function for [ss|ss] ERIs, the size is (1,1,1,1)
gabcd = 0;
%Screen if x for Interpolated_Boys will be large
g1 = basis_a.g(1);
g2 = basis_b.g(1);
g3 = basis_c.g(1);
g4 = basis_d.g(1);

% amax = g1.alpha;
% bmax = g2.alpha;
% cmax = g3.alpha;
% dmax = g4.alpha;
% pmax = amax+bmax;
% qmax = cmax+dmax;
% alphamax = pmax*qmax/(pmax+qmax);
% 
% Pxmax = (amax*g1.x0 + bmax*g2.x0)/pmax;
% Pymax = (amax*g1.y0 + bmax*g2.y0)/pmax;
% Pzmax = (amax*g1.z0 + bmax*g2.z0)/pmax;
% 
% Qxmax = (cmax*g3.x0 + dmax*g4.x0)/qmax;
% Qymax = (cmax*g3.y0 + dmax*g4.y0)/qmax;
% Qzmax = (cmax*g3.z0 + dmax*g4.z0)/qmax;
% 
% RPQ2max = (Pxmax-Qxmax)^2+(Pymax-Qymax)^2+(Pzmax-Qzmax)^2;


amin = basis_a.g(end).alpha;
bmin = basis_b.g(end).alpha;
cmin = basis_c.g(end).alpha;
dmin = basis_d.g(end).alpha;
pmin = amin+bmin;
qmin = cmin+dmin;
alphamin = pmin*qmin/(pmin+qmin);

Pxmin = (amin*g1.x0 + bmin*g2.x0)/pmin;
Pymin = (amin*g1.y0 + bmin*g2.y0)/pmin;
Pzmin = (amin*g1.z0 + bmin*g2.z0)/pmin;

Qxmin = (cmin*g3.x0 + dmin*g4.x0)/qmin;
Qymin = (cmin*g3.y0 + dmin*g4.y0)/qmin;
Qzmin = (cmin*g3.z0 + dmin*g4.z0)/qmin;

RPQ2min = (Pxmin-Qxmin)^2+(Pymin-Qymin)^2+(Pzmin-Qzmin)^2;

if (alphamin*RPQ2min > 20)

    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
%         g1 = basis_a.g(na);
%         aa = g1.alpha;
%         c1 = basis_a.c(na);
%         N1 = g1.N;
        
         na_data = ab_data{na,1};
         c1 = na_data.ca; 
         N1 = na_data.Na;
        
        temp = 0;
               
            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
%                 g2 = basis_b.g(nb);
%                 ab = g2.alpha;
%                 c2 = basis_b.c(nb);
%                 N2 = g2.N;
                
                nb_data = ab_data{na,nb};
                c2 = nb_data.cb; 
                N2 = nb_data.Nb;
                
%                 p = aa + ab;
%                 Px = (aa*g1.x0 + ab*g2.x0)/p;
%                 Py = (aa*g1.y0 + ab*g2.y0)/p;
%                 Pz = (aa*g1.z0 + ab*g2.z0)/p;
                
                p = nb_data.p;
                P = nb_data.P;
                
                %Not needed anymore
                %old RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %old RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                %RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

                
%                 rhoAB = aa*ab/p;
%                 Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                Kab = nb_data.Kab;
                
                temp2 = 0;
                
                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
%                         g3 = basis_c.g(nc);
%                         ac = g3.alpha;
%                         c3 = basis_c.c(nc);
%                         N3 = g3.N;
                        
                        %The nd index is still not defined, but
                        %pair_data{c,d} already contains the coefficients I
                        %need for the C function.
                        nc_data = cd_data{nc,1}; 
                        c3 = nc_data.ca; 
                        N3 = nc_data.Na;
                        
                        %old temp3 = zeros(Dim1,Dim2,Dim3,Dim4); %The same
                        %for all temp definitions
                        temp3 = 0;
                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
                            
%                             g4 = basis_d.g(nd);
%                             ad = g4.alpha;
%                             c4 = basis_d.c(nd);
%                             N4 = g4.N;
                            
                            nd_data = cd_data{nc,nd};
                            c4 = nd_data.cb;
                            N4 = nd_data.Nb;
                            
%                             q = ac+ad;
%                             Qx = (ac*g3.x0 + ad*g4.x0)/q;
%                             Qy = (ac*g3.y0 + ad*g4.y0)/q;
%                             Qz = (ac*g3.z0 + ad*g4.z0)/q;
                             q = nd_data.p;
                             Q = nd_data.P;

                            %Not needed anymore
                            %old RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %old RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
%                             RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
%                             
%                             
%                             
%                             rhoCD = ac*ad/q;
%                             Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
                            
                            Kcd = nd_data.Kab;
                            
                            
                            %RPQ = [P(1)-Qx;Py-Qy;Pz-Qz]; %column vector
                            RPQx = P(1)-Q(1);
                            RPQy = P(2)-Q(2);
                            RPQz = P(3)-Q(3);
                            %RPQ = P-Q;
                            RPQ2 = RPQx*RPQx+RPQy*RPQy+RPQz*RPQz;
                            alpha = q*p/(q+p);
                            
                                alphaRPQ2 = alpha*RPQ2;
                                Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));
                                %disp([num2str(na) num2str(nb) num2str(nc) num2str(nd) num2str(alphaRPQ2)]);
                                temp4 = Prefactor*0.886226925452758/sqrt(alphaRPQ2);

                            temp3 = temp3 + temp4*c4*N4;
                            
                        end
                        temp2 = temp2 + temp3*c3*N3;
                    end
                    temp = temp + temp2*c2*N2;
            end
            gabcd = gabcd + temp*c1*N1;
        
    end
    
else
    
    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
%         g1 = basis_a.g(na);
%         aa = g1.alpha;
%         c1 = basis_a.c(na);
%         N1 = g1.N;
        
         na_data = ab_data{na,1};
         c1 = na_data.ca; 
         N1 = na_data.Na;
        
        temp = 0;
               
            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
%                 g2 = basis_b.g(nb);
%                 ab = g2.alpha;
%                 c2 = basis_b.c(nb);
%                 N2 = g2.N;
                
                nb_data = ab_data{na,nb};
                c2 = nb_data.cb; 
                N2 = nb_data.Nb;
                
%                 p = aa + ab;
%                 Px = (aa*g1.x0 + ab*g2.x0)/p;
%                 Py = (aa*g1.y0 + ab*g2.y0)/p;
%                 Pz = (aa*g1.z0 + ab*g2.z0)/p;
                
                p = nb_data.p;
                P = nb_data.P;
                
                %Not needed anymore
                %old RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %old RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                %RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector

                
%                 rhoAB = aa*ab/p;
%                 Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                Kab = nb_data.Kab;
                
                temp2 = 0;
                
                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
%                         g3 = basis_c.g(nc);
%                         ac = g3.alpha;
%                         c3 = basis_c.c(nc);
%                         N3 = g3.N;
                        
                        %The nd index is still not defined, but
                        %pair_data{c,d} already contains the coefficients I
                        %need for the C function.
                        nc_data = cd_data{nc,1}; 
                        c3 = nc_data.ca; 
                        N3 = nc_data.Na;
                        
                        %old temp3 = zeros(Dim1,Dim2,Dim3,Dim4); %The same
                        %for all temp definitions
                        temp3 = 0;
                        for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
                            
%                             g4 = basis_d.g(nd);
%                             ad = g4.alpha;
%                             c4 = basis_d.c(nd);
%                             N4 = g4.N;
                            
                            nd_data = cd_data{nc,nd};
                            c4 = nd_data.cb;
                            N4 = nd_data.Nb;
                            
%                             q = ac+ad;
%                             Qx = (ac*g3.x0 + ad*g4.x0)/q;
%                             Qy = (ac*g3.y0 + ad*g4.y0)/q;
%                             Qz = (ac*g3.z0 + ad*g4.z0)/q;
                             q = nd_data.p;
                             Q = nd_data.P;

                            %Not needed anymore
                            %old RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %old RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
%                             RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
%                             
%                             
%                             
%                             rhoCD = ac*ad/q;
%                             Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
                            
                            Kcd = nd_data.Kab;
                            
                            
                            %RPQ = [P(1)-Qx;Py-Qy;Pz-Qz]; %column vector
                            RPQx = P(1)-Q(1);
                            RPQy = P(2)-Q(2);
                            RPQz = P(3)-Q(3);
                            %RPQ = P-Q;
                            RPQ2 = RPQx*RPQx+RPQy*RPQy+RPQz*RPQz;
                            alpha = q*p/(q+p);
                            
                            alphaRPQ2 = alpha*RPQ2;
                            Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));
    
                            tmp = Interpolated_Boys_2(alphaRPQ2,Boys_Table);


                            temp4 = Prefactor*tmp;
                            

                            temp3 = temp3 + temp4*c4*N4;
                            
                        end
                        temp2 = temp2 + temp3*c3*N3;
                    end
                    temp = temp + temp2*c2*N2;
            end
            gabcd = gabcd + temp*c1*N1;
        
    end    


%     for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
%         g1 = basis_a.g(na);
%         aa = g1.alpha;
%         c1 = basis_a.c(na);
%         N1 = g1.N;
%         
% %         na_data = ab_data{na,1};
% %         c1 = na_data.ca; 
% %         N1 = na_data.Na;
%         
%         temp = 0;
%                
%             for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
%                 g2 = basis_b.g(nb);
%                 ab = g2.alpha;
%                 c2 = basis_b.c(nb);
%                 N2 = g2.N;
%                 
% %                 nb_data = ab_data{na,nb};
% %                 c2 = nb_data.cb; 
% %                 N2 = nb_data.Nb;
%                 
%                 p = aa + ab;
%                 Px = (aa*g1.x0 + ab*g2.x0)/p;
%                 Py = (aa*g1.y0 + ab*g2.y0)/p;
%                 Pz = (aa*g1.z0 + ab*g2.z0)/p;
%                 
%                 %p = nb_data.p;
%                 %P = nb_data.P;
%                 
%                 %Not needed anymore
%                 %old RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
%                 %old RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
%                 RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector
% 
%                 
%                 rhoAB = aa*ab/p;
%                 Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
%                 %Kab = nb_data.Kab;
%                 
%                 temp2 = 0;
%                 
%                     for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function 
%                         g3 = basis_c.g(nc);
%                         ac = g3.alpha;
%                         c3 = basis_c.c(nc);
%                         N3 = g3.N;
%                         
%                         %The nd index is still not defined, but
%                         %pair_data{c,d} already contains the coefficients I
%                         %need for the C function.
% %                         nc_data = cd_data{nc,1}; 
% %                         c3 = nc_data.ca; 
% %                         N3 = nc_data.Na;
%                         
%                         %old temp3 = zeros(Dim1,Dim2,Dim3,Dim4); %The same
%                         %for all temp definitions
%                         temp3 = 0;
%                         for nd = 1:basis_d.n; %loops over number of primitives in 4th contracted basis functions.
%                             
%                             g4 = basis_d.g(nd);
%                             ad = g4.alpha;
%                             c4 = basis_d.c(nd);
%                             N4 = g4.N;
%                             
% %                             nd_data = cd_data{nc,nd};
% %                             c4 = nd_data.cb;
% %                             N4 = nd_data.Nb;
%                             
%                             q = ac+ad;
%                             Qx = (ac*g3.x0 + ad*g4.x0)/q;
%                             Qy = (ac*g3.y0 + ad*g4.y0)/q;
%                             Qz = (ac*g3.z0 + ad*g4.z0)/q;
% %                             q = nd_data.p;
% %                             Q = nd_data.P;
% 
%                             %Not needed anymore
%                             %old RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
%                             %old RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
%                             RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];
%                             
%                             
%                             
%                             rhoCD = ac*ad/q;
%                             Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));
%                             
%                             %Kcd = nd_data.Kab;
%                             
%                             RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
%                             %RPQ = P-Q;
%                             RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
%                             alpha = q*p/(q+p);
%                             
%                             alphaRPQ2 = alpha*RPQ2;
%                             Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(p+q));
%     
%                             %disp([num2str(na) ' ' num2str(nb) ' ' num2str(nc) ' ' num2str(nd) ' ' num2str(alphaRPQ2)]);
%                             try
%                             %tmp = Interpolated_Boys_0_small_x(alphaRPQ2,Boys_Table);
%                             tmp = Interpolated_Boys_2(alphaRPQ2,Boys_Table);
%                             catch
%                                 disp([num2str(aa) ' ' num2str(na) ' ' num2str(ab) ' ' num2str(nb) ' ' num2str(ac) ' ' num2str(nc) ' ' num2str(ad) ' ' num2str(nd) ' ' num2str(alphaRPQ2)]);
%                                 disp([num2str(amax) ' ' num2str(na) ' ' num2str(bmax) ' ' num2str(nb) ' ' num2str(cmax) ' ' num2str(nc) ' ' num2str(dmax) ' ' num2str(nd) ' ' num2str(alphamax*RPQ2max)]);
%                             end
%                             %tmp = Interpolated_Boys_0_small_x(alphaRPQ2,Boys_Table);
%                             temp4 = Prefactor*tmp;
%                             
%                             temp3 = temp3 + temp4*c4*N4;
%                             
%                         end
%                         temp2 = temp2 + temp3*c3*N3;
%                     end
%                     temp = temp + temp2*c2*N2;
%             end
%             gabcd = gabcd + temp*c1*N1;
%         
%     end
% 
% end


end