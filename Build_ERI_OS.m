function gabcd = Build_ERI_OS(basis,Shell_List,Boys_Table,pair_data)
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
                            gSSSSNValues = primitiveFactorsSSSS_2(basis_a,basis_b,basis_c,basis_d,Boys_Table);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = sum(gSSSSNValues(:,1));

                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_2(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_vrr(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
                            %permuted_LSSS = contract_vrr(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz); %b and a are permuted
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [2 1 3 4]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
                            permuted_LSSS = contract_vrr(basis_c,basis_d,basis_a,basis_b,basis_c.L,basis_d.L,basis_a.L,basis_b.L,Boys_Table,nz);
                            permuted_LSSS = permute(permuted_LSSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [1 2 3 4]);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
                            permuted_LSSS = contract_vrr(basis_d,basis_c,basis_a,basis_b,basis_d.L,basis_c.L,basis_a.L,basis_b.L,Boys_Table,nz); %b and a are permuted
                            permuted_LSSS = permute(permuted_LSSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [1 2 4 3]);


                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);


                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        else
                            break;
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);
                        end

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        gabcd(intnu,intmu,intka,intla) = permute(gabcd(intmu,intnu,intka,intla),[2 1 3 4]);
                        gabcd(intmu,intnu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[1 2 4 3]);
                        gabcd(intnu,intmu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[2 1 4 3]);
                        gabcd(intka,intla,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 1 2]);
                        gabcd(intla,intka,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 1 2]);
                        gabcd(intka,intla,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 2 1]);
                        gabcd(intla,intka,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 2 1]);

                end
            else
                for d = 1:c
                    la_begin = Shell_List(d,1);
                    la_end = Shell_List(d,2);
                    %ad = Shell_List(d,2);
                    basis_d = basis{d};
                    cd_data = pair_data{c,d};

                        if (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %SS|SS integrals
                            gSSSSNValues = primitiveFactorsSSSS_2(basis_a,basis_b,basis_c,basis_d,Boys_Table);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = sum(gSSSSNValues(:,1));

                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_2(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_ssss_3(basis_a,basis_b,basis_c,basis_d,Boys_Table,ab_data,cd_data);
                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L == 0) %LS|SS integrals
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = contract_vrr(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table,nz);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);
                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %SL|SS integrals
                            %permuted_LSSS = contract_vrr(basis_b,basis_a,basis_c,basis_d,basis_b.L,basis_a.L,basis_c.L,basis_d.L,Boys_Table,nz); %b and a are permuted
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [2 1 3 4]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %SS|LS integrals
                            permuted_LSSS = contract_vrr(basis_c,basis_d,basis_a,basis_b,basis_c.L,basis_d.L,basis_a.L,basis_b.L,Boys_Table,nz);
                            permuted_LSSS = permute(permuted_LSSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [1 2 3 4]);
                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %SS|SL integrals
                            permuted_LSSS = contract_vrr(basis_d,basis_c,basis_a,basis_b,basis_d.L,basis_c.L,basis_a.L,basis_b.L,Boys_Table,nz); %b and a are permuted
                            permuted_LSSS = permute(permuted_LSSS, [3 4 1 2]);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = permute(permuted_LSSS, [1 2 4 3]);


                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L == 0) %LaS|LcS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %SLb|LcS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L == 0 && basis_d.L > 0) %LaS|SLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %SLb|SLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L == 0) %LaLb|SS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);


                        elseif (basis_a.L == 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %SS|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L == 0) %LaLb|LcS integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L == 0 && basis_d.L > 0) %LaLb|SLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L == 0 && basis_c.L > 0 && basis_d.L > 0) %LaS|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L == 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %SLb|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        elseif (basis_a.L > 0 && basis_b.L > 0 && basis_c.L > 0 && basis_d.L > 0) %LaLb|LcLd integrals

                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = shellOS(basis_a,basis_b,basis_c,basis_d,basis_a.L,basis_b.L,basis_c.L,basis_d.L,Boys_Table);

                        else
                            break;
                            %gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);
                        end

                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        gabcd(intnu,intmu,intka,intla) = permute(gabcd(intmu,intnu,intka,intla),[2 1 3 4]);
                        gabcd(intmu,intnu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[1 2 4 3]);
                        gabcd(intnu,intmu,intla,intka) = permute(gabcd(intmu,intnu,intka,intla),[2 1 4 3]);
                        gabcd(intka,intla,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 1 2]);
                        gabcd(intla,intka,intmu,intnu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 1 2]);
                        gabcd(intka,intla,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[3 4 2 1]);
                        gabcd(intla,intka,intnu,intmu) = permute(gabcd(intmu,intnu,intka,intla),[4 3 2 1]);

                end
            end
        end
    end
end

% for intmu = 1:Ncont
%     for intnu = 1:intmu
%         for intka = 1:intmu
%            if (intka == intmu)
%                 for intla = 1:intnu
%                         gabcd(intnu,intmu,intka,intla) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intmu,intnu,intla,intka) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intnu,intmu,intla,intka) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intka,intla,intmu,intnu) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intla,intka,intmu,intnu) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intka,intla,intnu,intmu) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intla,intka,intnu,intmu) = gabcd(intmu,intnu,intka,intla);
%                 end
%            else
%                 for intla = 1:intka
%                         gabcd(intnu,intmu,intka,intla) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intmu,intnu,intla,intka) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intnu,intmu,intla,intka) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intka,intla,intmu,intnu) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intla,intka,intmu,intnu) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intka,intla,intnu,intmu) = gabcd(intmu,intnu,intka,intla);
%                         gabcd(intla,intka,intnu,intmu) = gabcd(intmu,intnu,intka,intla);
%                 end
%            end
%         end
%     end
% end



end
