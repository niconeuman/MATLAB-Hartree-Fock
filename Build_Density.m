function D = Build_Density(eigenvecs,Nel)

nb = size(eigenvecs,1);
D = zeros(nb,nb);

% for n = 1:nb
%     for m = 1:nb
%         for i = 1:Nel/2
%             D(n,m) = D(n,m)+2*eigenvecs(n,i)*eigenvecs(m,i);
%         end
%     end
% end

for i = 1:Nel/2
    D = D + 2*eigenvecs(:,i)*eigenvecs(:,i)';
end