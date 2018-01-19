function gabcd = Build_gabcd_sym(basis,Shell_List)

Ncont = Shell_List(end,2);
nb = size(basis,1);
gabcd = zeros(Ncont,Ncont,Ncont,Ncont);

for a = 1:nb
    mu_begin = Shell_List(a,1);
    mu_end = Shell_List(a,2);
    %aa = Shells_List(a,2);
    basis_a = basis{a};
    Dim1 = (basis_a.L+1)*(basis_a.L+2)/2;
    for b = 1:a
        nu_begin = Shell_List(b,1);
        nu_end = Shell_List(b,2);
        %ab = Shells_List(b,2);
        basis_b = basis{b};
        Dim2 = (basis_b.L+1)*(basis_b.L+2)/2;
        for c = 1:a
        ka_begin = Shell_List(c,1);
        ka_end = Shell_List(c,2);
        %ac = Shells_List(c,2);
        basis_c = basis{c};
        Dim3 = (basis_c.L+1)*(basis_c.L+2)/2;
            if (c == a)
                for d = 1:b
                    la_begin = Shell_List(d,1);
                    la_end = Shell_List(d,2);
                    %ad = Shells_List(d,2);
                    basis_d = basis{d};
                    Dim4 = (basis_d.L+1)*(basis_d.L+2)/2;    
                        
                    
                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        
                        gabcd(intmu,intnu,intka,intla) = gabcd(intmu,intnu,intka,intla) + ones(Dim1,Dim2,Dim3,Dim4);
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
                    %ad = Shells_List(d,2);
                    basis_d = basis{d};
                    Dim4 = (basis_d.L+1)*(basis_d.L+2)/2; 
                    
                        intmu = mu_begin:mu_end;
                        intnu = nu_begin:nu_end;
                        intka = ka_begin:ka_end;
                        intla = la_begin:la_end;
                        gabcd(intmu,intnu,intka,intla) = gabcd(intmu,intnu,intka,intla) + ones(Dim1,Dim2,Dim3,Dim4);
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


end