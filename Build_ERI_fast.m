function [gabcd,int_error] = Build_ERI_fast(basis,Shell_List,Boys_Table,pair_data,S,T,Pcell,pcell,gabcd_exact)

Ncont = Shell_List(end,2);
nb = size(basis,1);
gabcd = zeros(Ncont,Ncont,Ncont,Ncont);
int_error = zeros(Ncont,Ncont,Ncont,Ncont);

%This will give me the nonzero elements that contribute in the [a-2s|ss]
%terms in the vertical recursion relations.
nz = cell(1,10);
nz{2} = [1;0;0;1;0;1];
nz{3} = [2;1;1;0;0;0;2;1;0;2];
nz{4} = [3;2;2;1;1;1;0;0;0;0;3;2;1;0;3];
nz{5} = [4;3;3;2;2;2;1;1;1;1;0;0;0;0;0;4;3;2;1;0;4];
nz{6} = [5;4;4;3;3;3;2;2;2;2;1;1;1;1;1;0;0;0;0;0;0;5;4;3;2;1;0;5];
nz{7} = [6;5;5;4;4;4;3;3;3;3;2;2;2;2;2;1;1;1;1;1;1;0;0;0;0;0;0;0;6;5;4;3;2;1;0;6];
nz{8} = [7;6;6;5;5;5;4;4;4;4;3;3;3;3;3;2;2;2;2;2;2;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;7;6;5;4;3;2;1;0;7];
nz{9} = [8;7;7;6;6;6;5;5;5;5;4;4;4;4;4;3;3;3;3;3;3;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;8;7;6;5;4;3;2;1;0;8];
nz{10} = [9;8;8;7;7;7;6;6;6;6;5;5;5;5;5;4;4;4;4;4;4;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;9;8;7;6;5;4;3;2;1;0;9];

unique = [2, 4,5, 7,8,9, 11,12,13,14, 16,17,18,19,20, 22,23,24,25,26,27,  29,30,31,32,33,34,35, 37,38,39,40,41,42,43,44, 46,47,48,49,50,51,52,53,54, 56,57,58,59,60,61,62,63,64,65 ]; %up to [ms|ss]

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
                            
                            
                            alpha = (pcell{a,b}*pcell{c,d})/(pcell{a,b}+pcell{c,d});
                            x = alpha*sum((Pcell{a,b}-Pcell{c,d}).^2);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = 2*sqrt(alpha)/sqrt(pi)*S(a,b)*S(c,d)*Interpolated_Boys_2(x,Boys_Table);
                            int_error(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd_exact(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end)-gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);
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
                            
                            
                            alpha = (pcell{a,b}*pcell{c,d})/(pcell{a,b}+pcell{c,d});
                            x = alpha*sum((Pcell{a,b}-Pcell{c,d}).^2);
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = 2*sqrt(alpha)/sqrt(pi)*S(a,b)*S(c,d)*Interpolated_Boys_2(x,Boys_Table);
                            int_error(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd_exact(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end)-gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);
                        
                        else
                            gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end) = gabcd(mu_begin:mu_end,nu_begin:nu_end,ka_begin:ka_end,la_begin:la_end);


                        end
                end
            end
        end
    end
end

for a = 1:nb
    for b = 1:nb
        for c = 1:nb
            for d = 1:nb
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