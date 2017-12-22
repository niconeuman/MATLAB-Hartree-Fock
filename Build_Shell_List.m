function Shell_List = Build_Shell_List(basis)
%This function builds a shell list which has 3 columns, 
%[mu_begin mu_end nb]
%It is a much smaller version of Build_Shells, but contains the same
%information in a non-redundant way.

nb = size(basis,1);
Shell_List = zeros(nb,3);
mu_begin = 1;

for t = 1:nb
    Length_mu = (2*basis{t}.L+1);
    mu_end = mu_begin + Length_mu - 1;
    Shell_List(t,:) = [mu_begin mu_end t];
    mu_begin = mu_end + 1;
end


end