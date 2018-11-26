function gabcd = contract_ssss_2(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data)

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

%As it is a function for [ss|ss] ERIs, the size is (1,1,1,1);
gabcd = 0;



    for na=1:basis_a.n %loops over the number of primitives in the 1st contracted basis function
        %old g1 = basis_a.g(na);
        %old aa = g1.alpha;
        %old c1 = basis_a.c(na);
        %old N1 = g1.N;

        na_data = ab_data{na,1};
        c1 = na_data.ca;
        N1 = na_data.Na;

        temp = 0;

            for nb=1:basis_b.n %loops over number of primitives in the 2nd contracted basis function
                %old g2 = basis_b.g(nb);
                %old ab = g2.alpha;
                %old c2 = basis_b.c(nb);
                %old N2 = g2.N;

                nb_data = ab_data{na,nb};
                c2 = nb_data.cb;
                N2 = nb_data.Nb;

                %old p = aa + ab;
                %old Px = (aa*g1.x0 + ab*g2.x0)/p;
                %old Py = (aa*g1.y0 + ab*g2.y0)/p;
                %old Pz = (aa*g1.z0 + ab*g2.z0)/p;

                p = nb_data.p;
                P = nb_data.P;

                %Not needed anymore
                %old RPA = [Px-g1.x0;Py-g1.y0;Pz-g1.z0];
                %old RPB = [Px-g2.x0;Py-g2.y0;Pz-g2.z0];
                %old RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0]; %column vector


                %old rhoAB = aa*ab/p;
                %old Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                Kab = nb_data.Kab;

                temp2 = 0;

                    for nc = 1:basis_c.n; %loop over number of primitives in 3rd contracted basis function
                        %old g3 = basis_c.g(nc);
                        %old c = g3.alpha;
                        %old c3 = basis_c.c(nc);
                        %old N3 = g3.N;

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

                            %old g4 = basis_d.g(nd);
                            %old ad = g4.alpha;
                            %old c4 = basis_d.c(nd);
                            %old N4 = g4.N;

                            nd_data = cd_data{nc,nd};
                            c4 = nd_data.cb;
                            N4 = nd_data.Nb;

                            %old q = ac+ad;
                            %old Qx = (ac*g3.x0 + ad*g4.x0)/q;
                            %old Qy = (ac*g3.y0 + ad*g4.y0)/q;
                            %old Qz = (ac*g3.z0 + ad*g4.z0)/q;
                            q = nd_data.p;
                            Q = nd_data.P;

                            %Not needed anymore
                            %old RQC = [Qx-g3.x0;Qy-g3.y0;Qz-g3.z0];
                            %old RQD = [Qx-g4.x0;Qy-g4.y0;Qz-g4.z0];
                            %old RCD = [g3.x0-g4.x0;g3.y0-g4.y0;g3.z0-g4.z0];



                            %old rhoCD = ac*ad/q;
                            %old Kcd = exp(-rhoCD*(RCD(1)^2+RCD(2)^2+RCD(3)^2));

                            Kcd = nd_data.Kab;

                            %old RPQ = [Px-Qx;Py-Qy;Pz-Qz]; %column vector
                            RPQ = [P(1)-Q(1);P(2)-Q(2);P(3)-Q(3)];
                            RPQ2 = RPQ(1)^2+RPQ(2)^2+RPQ(3)^2;
                            alpha = q*p/(q+p);



                            temp3 = temp3 + g_SSSS_0(Kab,Kcd,p,q,alpha,RPQ2,Boys_Table)*c4*N4;

                        end
                        temp2 = temp2 + temp3*c3*N3;
                    end
                    temp = temp + temp2*c2*N2;
            end
            gabcd = gabcd + temp*c1*N1;

    end




end
