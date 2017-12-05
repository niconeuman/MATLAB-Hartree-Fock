function [NonDiag,Ratio,Schwarz_bound,Diag_ab,Diag_cd] = plot_gabcd(gabcd)

Ncont = size(gabcd,1);

t = 1;
for a = 1:Ncont
    for b = 1:a
        for c = 1:Ncont
            for d = 1:c
                NonDiag(t) = gabcd(a,b,c,d);
                Diag_ab(t) = gabcd(a,b,a,b);
                Diag_cd(t) = gabcd(c,d,c,d);
                Schwarz_bound(t) = sqrt(gabcd(a,b,a,b))*sqrt(gabcd(c,d,c,d));
                Ratio(t) = gabcd(a,b,c,d)/Schwarz_bound(t);
                t = t+1;
            end
        end
    end
end



end