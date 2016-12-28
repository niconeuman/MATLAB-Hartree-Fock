function [Shells,NShells] = Build_Shells(basis)
%This function creates shells which are arrays 
nb = size(basis,2);

%Something is wrong with this function

%Structure of the Shells matrix = [mu_begin mu_end a nu_begin nu_end b kappa_begin kappa_end c lambda_begin lambda_end d] 
Shells = zeros(2000,12); %2000 shells is more than enough for now, but maybe later I have to change it
t = 1;
mu_begin = 1;
nu_begin = 1;
kappa_begin = 1;
lambda_begin = 1;

%There is a subtlety. (Contracted) Basis functions are either s, p, d, f.
%They don't specify each of the px, py, pz, components.
%But the gabcd matrix and other matrices contain distinct elements for the
%<px s|s s>, <py s|s s>, <pz s|s s>, etc. So the mu_begin and mu_end
%indices, etc, will not coincide with the a, b, c, d used for calling basis
%functions. Have to be careful with this.

for a = 1:nb
    for b = 1:nb
        for c = 1:nb
            for d = 1:nb
                Length_mu = (2*basis{a}.L+1); %if L = 0,1,2, etc Length_mu = 1,3,5, etc
                Length_nu = (2*basis{b}.L+1);  
                Length_kappa = (2*basis{c}.L+1);  
                Length_lambda = (2*basis{d}.L+1);
                mu_end = mu_begin + Length_mu - 1;
                nu_end = nu_begin + Length_nu - 1;
                kappa_end = kappa_begin + Length_kappa - 1;
                lambda_end = lambda_begin + Length_lambda - 1;
                
                
                Shells(t,:) = [mu_begin mu_end a nu_begin nu_end b kappa_begin kappa_end c lambda_begin lambda_end d];
                t = t + 1;
                mu_begin = mu_end + 1;
                nu_begin = nu_end + 1;
                kappa_begin = kappa_end + 1;
                lambda_begin = lambda_end + 1;
            end
        end
    end
end
Shells = Shells(1:(t-1),:);
NShells = size(Shells,1);
end