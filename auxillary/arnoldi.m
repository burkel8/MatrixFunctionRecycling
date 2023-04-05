%Function to apply Arnoldi algorithm to build basis for Krylov subspace Km(A,b)

%Input: Aop: function operator to apply appropriate matrix to a vector
%         b: vector b from f(A)b
%         p: problem struct
%Output:  H,V - Hessenberg matrix and Arnoldi basis built from Arnoldi process. 

function [V,H] = arnoldi(Aop, b, p)

V = zeros(p.n,p.m+1);
H = zeros(p.m+1,p.m);

V(:,1) = b/norm(b);

for j=1:p.m

   V(:,j+1) = Aop(V(:,j));
  
   %modified GS
   for i=1:j
       H(i,j)= V(:,i)'*V(:,j+1);
       V(:,j+1)= V(:,j+1) - V(:,i)*H(i,j);
   end
    
    H(j+1,j) = norm(V(:,j+1));
    V(:,j+1) = V(:,j+1)/H(j+1,j);
end

end
