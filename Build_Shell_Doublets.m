function [Shell_Doublets,NShell_Doublets] = Build_Shell_Doublets(basis)

%Something was wrong. mu needs to loop over a, not a,b. nu loops over b.
%When I increase mu, I need to reset nu
nb = size(basis,1);

%Structure of the Shells matrix = [mu_begin mu_end a nu_begin nu_end b ] 
Shell_Doublets = zeros(2000,6); %2000 shells is more than enough for now, but maybe later I have to change it
t = 1;
mu_begin = 1;
nu_begin = 1;

for a = 1:nb
    Length_mu = (2*basis{a}.L+1); %if L = 0,1,2, etc Length_mu = 1,3,5, etc
    mu_end = mu_begin + Length_mu - 1;
    for b = 1:nb  
        Length_nu = (2*basis{b}.L+1);
        nu_end = nu_begin + Length_nu - 1; 
        Shell_Doublets(t,:) = [mu_begin mu_end a nu_begin nu_end b ];
        t = t + 1;
        
        nu_begin = nu_end + 1;
        
    end
    mu_begin = mu_end + 1;
    nu_begin = 1;
end
Shell_Doublets = Shell_Doublets(1:(t-1),:);
NShell_Doublets = size(Shell_Doublets,1);

end