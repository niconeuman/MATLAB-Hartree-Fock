function gabcd = Build_ERI_2(basis,Shell_List,Boys_Table,pair_data)
%2 jul 2017. This is a function similar to Build_ERI, but which uses 4
%loops because it uses the data of Shell_List (an nb x 3 matrix), instead
%of Shells, (an nb^4 x 12) matrix.

%Size of the contracted basis list (size of gabcd matrix)
Ncont = Shell_List(end,2);
nb = size(basis,1);
gabcd = zeros(Ncont,Ncont,Ncont,Ncont);

%This will give me the nonzero elements that contribute in the [a-2s|ss]
%terms in the vertical recursion relations.
nz = cell(1,12);
nz{2} = [1;0;0;1;0;1];
nz{3} = [2;1;1;0;0;0;2;1;0;2];
nz{4} = [3;2;2;1;1;1;0;0;0;0;3;2;1;0;3];
nz{5} = [4;3;3;2;2;2;1;1;1;1;0;0;0;0;0;4;3;2;1;0;4];
nz{6} = [5;4;4;3;3;3;2;2;2;2;1;1;1;1;1;0;0;0;0;0;0;5;4;3;2;1;0;5];
nz{7} = [6;5;5;4;4;4;3;3;3;3;2;2;2;2;2;1;1;1;1;1;1;0;0;0;0;0;0;0;6;5;4;3;2;1;0;6];
nz{8} = [7;6;6;5;5;5;4;4;4;4;3;3;3;3;3;2;2;2;2;2;2;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;7;6;5;4;3;2;1;0;7];
nz{9} = [8;7;7;6;6;6;5;5;5;5;4;4;4;4;4;3;3;3;3;3;3;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;8;7;6;5;4;3;2;1;0;8];
nz{10} = [9;8;8;7;7;7;6;6;6;6;5;5;5;5;5;4;4;4;4;4;4;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;9;8;7;6;5;4;3;2;1;0;9];
nz{11} = [10;9;9;8;8;8;7;7;7;7;6;6;6;6;6;5;5;5;5;5;5;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;10;9;8;7;6;5;4;3;2;1;0;10];
nz{12} = [11;10;10;9;9;9;8;8;8;8;7;7;7;7;7;6;6;6;6;6;6;5;5;5;5;5;5;5;4;4;4;4;4;4;4;4;3;3;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;11;10;9;8;7;6;5;4;3;2;1;0;11];

unique = [2, 4,5, 7,8,9, 11,12,13,14, 16,17,18,19,20, 22,23,24,25,26,27,  29,30,31,32,33,34,35, 37,38,39,40,41,42,43,44, 46,47,48,49,50,51,52,53,54, 56,57,58,59,60,61,62,63,64,65, 67,68,69,70,71,72,73,74,75,76,77, 79,80,81,82,83,84,85,86,87,88,89,90]; %up to [ms|ss]

for a = 1:nb
    mu_begin = Shell_List(a,1);
    mu_end = Shell_List(a,2);
    %aa = Shells_List(a,2);
    basis_a = basis{a};
    for b = 1:a
        nu_begin = Shell_List(b,1);
        nu_end = Shell_List(b,2);
        %ab = Shells_List(b,2);
        basis_b = basis{b};
        ab_data = pair_data{a,b};
        for c = 1:a
        ka_begin = Shell_List(c,1);
        ka_end = Shell_List(c,2);
        %ac = Shells_List(c,2);
        basis_c = basis{c};
            if (c == a)
                for d = 1:b
                    la_begin = Shell_List(d,1);
                    la_end = Shell_List(d,2);
                    %ad = Shells_List(d,2);
                    basis_d = basis{d};
                    cd_data = pair_data{c,d};
                
                        if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_2(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_vrr(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz);
     
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
                            permuted_LSSS = contract_vrr(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz); %b and a are permuted
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [2 1 3 4]);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
                            permuted_LSSS = contract_vrr(basis_c,basis_b,basis_a,basis_d,basis_c.L,basis_b.L,basis_a.L,basis_d.L,Boys_Table,nz); %b and a are permuted
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [3 2 1 4]);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
                            permuted_LSSS = contract_vrr(basis_d,basis_b,basis_c,basis_a,basis_d.L,basis_b.L,basis_c.L,basis_a.L,Boys_Table,nz); %b and a are permuted
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [4 2 3 1]);
    
         
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_c.L > basis_a.L)
                                LaSLcS = contract_ascs(basis_c,basis_b,basis_a,basis_d,basis_c.L,basis_b.L,basis_a.L,basis_d.L,Boys_Table,nz,unique);
                                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                %This would not be OK,
                                %check!!!!!!!!!!!!!!!!!!!!!!
                                LaSLcS = permute(LaSLcS, [3 2 1 4]);
                            else
                                LaSLcS = contract_ascs(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            
                            %LaSLcS is structured as a [DimA, DimC,1,1] matrix, for convenience in matrix operations.
                            LaSLcS = permute(LaSLcS, [1 3 2 4]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_c.L > basis_b.L)
                                LaSLcS = contract_ascs(basis_c,basis_a,basis_b,basis_d,basis_c.L,basis_a.L,basis_b.L,basis_d.L,Boys_Table,nz,unique);
                                LaSLcS = permute(LaSLcS, [3 2 1 4]);
                            else
                                LaSLcS = contract_ascs(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            
                            %LaSLcS is structured as a [DimA, DimC,1,1] matrix, for convenience in matrix operations.
                            LaSLcS = permute(LaSLcS, [3 1 2 4]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;                            
                            
                            
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_d.L > basis_a.L)
                                LaSLcS = contract_ascs(basis_d,basis_b,basis_a,basis_c,basis_d.L,basis_b.L,basis_a.L,basis_c.L,Boys_Table,nz,unique);
                                LaSLcS = permute(LaSLcS, [3 2 1 4]);
                            else
                                LaSLcS = contract_ascs(basis_a,basis_b,basis_d,basis_c,basis_a.L,basis_b.L,basis_d.L,basis_c.L,Boys_Table,nz,unique);
                            end
                            
                            %[DimA,DimD,1,1] -> [DimA,1,1,DimD]
                            LaSLcS = permute(LaSLcS, [1 4 3 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;
                            
                                %     elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals
                                %I need the hrr
                                %     elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals
                                %I need the hrr and permute
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals    
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_d.L > basis_b.L)
                                LaSLcS = contract_ascs(basis_d,basis_a,basis_b,basis_c,basis_d.L,basis_a.L,basis_b.L,basis_c.L,Boys_Table,nz,unique);
                                LaSLcS = permute(LaSLcS, [3 2 1 4]); %This is always so
                            else
                                LaSLcS = contract_ascs(basis_b,basis_a,basis_d,basis_c,basis_b.L,basis_a.L,basis_d.L,basis_c.L,Boys_Table,nz,unique);
                            end
                            
                            %This changes every time
                            %LaSLcS is structured as a [DimB,1,DimD,1]
                            LaSLcS = permute(LaSLcS, [2 1 4 3]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;                            
                        
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals   
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                            
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals   
                            if basis_d.L > basis_c.L
                                LaLbSS = hrr2D(basis_d,basis_c,basis_a,basis_b,basis_d.L,basis_c.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_c,basis_d,basis_a,basis_b,basis_c.L,basis_d.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                            end
                                LaLbSS = permute(LaLbSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                            
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;                            
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals
                            if basis_d.L > basis_c.L
                                LaLbSS = hrr2D(basis_d,basis_c,basis_a,basis_b,basis_d.L,basis_c.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_c,basis_d,basis_a,basis_b,basis_c.L,basis_d.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                            end
                                LaLbSS = permute(LaLbSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;  
                        else
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);


                        end

    
    
                
                
                
                end
            else
                for d = 1:c
                    la_begin = Shell_List(d,1);
                    la_end = Shell_List(d,2);
                    %ad = Shell_List(d,2);
                    basis_d = basis{d};
                    cd_data = pair_data{c,d};
                    
                        if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_2(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_vrr(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz);
     
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
                            permuted_LSSS = contract_vrr(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz); %b and a are permuted
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [2 1 3 4]);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
                            permuted_LSSS = contract_vrr(basis_c,basis_b,basis_a,basis_d,basis_c.L,basis_b.L,basis_a.L,basis_d.L,Boys_Table,nz); %b and a are permuted
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [3 2 1 4]);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
                            permuted_LSSS = contract_vrr(basis_d,basis_b,basis_c,basis_a,basis_d.L,basis_b.L,basis_c.L,basis_a.L,Boys_Table,nz); %b and a are permuted
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [4 2 3 1]);
    
         
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_c.L > basis_a.L)
                                LaSLcS = contract_ascs(basis_c,basis_b,basis_a,basis_d,basis_c.L,basis_b.L,basis_a.L,basis_d.L,Boys_Table,nz,unique);
                                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                %This would not be OK,
                                %check!!!!!!!!!!!!!!!!!!!!!!
                                LaSLcS = permute(LaSLcS, [3 2 1 4]);
                            else
                                LaSLcS = contract_ascs(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            
                            %LaSLcS is structured as a [DimA, DimC,1,1] matrix, for convenience in matrix operations.
                            LaSLcS = permute(LaSLcS, [1 3 2 4]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_c.L > basis_b.L)
                                LaSLcS = contract_ascs(basis_c,basis_a,basis_b,basis_d,basis_c.L,basis_a.L,basis_b.L,basis_d.L,Boys_Table,nz,unique);
                                LaSLcS = permute(LaSLcS, [3 2 1 4]);
                            else
                                LaSLcS = contract_ascs(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            
                            %LaSLcS is structured as a [DimA, DimC,1,1] matrix, for convenience in matrix operations.
                            LaSLcS = permute(LaSLcS, [3 1 2 4]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;                            
                            
                            
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_d.L > basis_a.L)
                                LaSLcS = contract_ascs(basis_d,basis_b,basis_a,basis_c,basis_d.L,basis_b.L,basis_a.L,basis_c.L,Boys_Table,nz,unique);
                                LaSLcS = permute(LaSLcS, [3 2 1 4]);
                            else
                                LaSLcS = contract_ascs(basis_a,basis_b,basis_d,basis_c,basis_a.L,basis_b.L,basis_d.L,basis_c.L,Boys_Table,nz,unique);
                            end
                            
                            %[DimA,DimD,1,1] -> [DimA,1,1,DimD]
                            LaSLcS = permute(LaSLcS, [1 4 3 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;
                            
                                %     elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals
                                %I need the hrr


                        
                        
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals    
                            %The following conditional should take care of cases
                            %where Lc > La, correctly
                            if (basis_d.L > basis_b.L)
                                LaSLcS = contract_ascs(basis_d,basis_a,basis_b,basis_c,basis_d.L,basis_a.L,basis_b.L,basis_c.L,Boys_Table,nz,unique);
                                LaSLcS = permute(LaSLcS, [3 2 1 4]); %This is always so
                            else
                                LaSLcS = contract_ascs(basis_b,basis_a,basis_d,basis_c,basis_b.L,basis_a.L,basis_d.L,basis_c.L,Boys_Table,nz,unique);
                            end
                            
                            %This changes every time
                            %LaSLcS is structured as a [DimB,1,DimD,1]
                            LaSLcS = permute(LaSLcS, [2 1 4 3]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaSLcS;                            
                        
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals   
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals   
                            if basis_d.L > basis_c.L
                                LaLbSS = hrr2D(basis_d,basis_c,basis_a,basis_b,basis_d.L,basis_c.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_c,basis_d,basis_a,basis_b,basis_c.L,basis_d.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                            end
                                LaLbSS = permute(LaLbSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;    
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals
                            if basis_b.L > basis_a.L
                                LaLbSS = hrr2D(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz,unique);
                            end
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;                            
                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals
                            if basis_d.L > basis_c.L
                                LaLbSS = hrr2D(basis_d,basis_c,basis_a,basis_b,basis_d.L,basis_c.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                                LaLbSS = permute(LaLbSS, [2 1 3 4]);
                            else
                                LaLbSS = hrr2D(basis_c,basis_d,basis_a,basis_b,basis_c.L,basis_d.L,basis_a.L,basis_b.L,Boys_Table,nz,unique);
                            end
                                LaLbSS = permute(LaLbSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = LaLbSS;  
                        else
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);


                        end

                end
            end
        end
    end
end

for a = 1:Ncont
    for b = 1:Ncont
        for c = 1:Ncont
            for d = 1:Ncont
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(b,a,c,d);
                end
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(a,b,d,c);
                end
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(b,a,d,c);
                end
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(d,c,a,b);
                end
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(d,c,b,a);
                end
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(c,d,a,b);
                end
                if (gabcd(a,b,c,d) == 0)
                    gabcd(a,b,c,d) = gabcd(c,d,b,a);
                end
            end
        end
    end
end


end