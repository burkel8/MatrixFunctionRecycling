function [U,D] = scale_cols_of_U(U,k)
D = zeros(k,k);
for i = 1:k
   D(i,i) = 1.0/norm(U(:,i));
   U(:,i) = U(:,i)*D(i,i);
end
end