function pair_data = basis_products(basis,Shell_Doublets)

%This function calculates shell pair information for further use in the two
%electron integral loops.

%The function takes a very small amount of time compared with the electron
%repulsion and even the one electron integrals, so there is not much point
%in optimizing it.
Ncont = Shell_Doublets(end,2); %This is the final mu_end, number of contracted basis functions

pair_data = cell(Ncont,Ncont); %I need a cell structure because I have many fields

nb = size(basis,1);

for a = 1:nb %loop over number of contracted functions
    for b = 1:nb
        pair_data{a,b} = cell(basis{a}.n,basis{a}.n);
        for nba = 1:basis{a}.n
            aa = basis{a}.g(nba).alpha;
            A = [basis{a}.g(nba).x0;basis{a}.g(nba).y0;basis{a}.g(nba).z0];
            Na = basis{a}.g(nba).N;
            ca = basis{a}.c(nba);
            for nbb = 1:basis{b}.n
                ab = basis{b}.g(nbb).alpha;
                p = aa+ab;
                B = [basis{b}.g(nbb).x0;basis{b}.g(nbb).y0;basis{b}.g(nbb).z0];
                Nb = basis{b}.g(nbb).N;
                cb = basis{b}.c(nbb);
                P = (aa*A+ab*B)/p; 
                RAB = B-A;
                RPA = P-A;
                RPB = P-B;
                rhoAB = aa*ab/p;
                Kab = exp(-rhoAB*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                
                    
                    pair_data{a,b}{nba,nbb}.RAB = RAB;
                    pair_data{a,b}{nba,nbb}.P = P;
                    pair_data{a,b}{nba,nbb}.p = p;
                    pair_data{a,b}{nba,nbb}.RPA = RPA;
                    pair_data{a,b}{nba,nbb}.RPB = RPB;
                    pair_data{a,b}{nba,nbb}.rhoAB = rhoAB;
                    pair_data{a,b}{nba,nbb}.Kab = Kab;
                    pair_data{a,b}{nba,nbb}.Na = Na;
                    pair_data{a,b}{nba,nbb}.Nb = Nb;
                    pair_data{a,b}{nba,nbb}.ca = ca;
                    pair_data{a,b}{nba,nbb}.cb = cb;
            end
        end
    end
end

end