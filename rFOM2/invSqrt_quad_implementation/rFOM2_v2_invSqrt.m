%Function to compute approximation to f(A)b when f is the inverse square root function
% using implementation 2 of r(FOM)^2. The quadrature rule used is outlined
% in section 7 of the preprint.

%Inputs:
% A struct p with the following attributes
%         b - vector b for which we want to approximate f(A)b
%        V,H - basis of Krylov subspace and hessenberg matrix H constructed
%        from Arnoldi process to build basis for Km(A,b)
%        m - dimension of Krylov subspace
%        k - dimension of recycling subspace
%        U - recycling subspace
%        C - matrix C = A*U.
%        num_quad - number of quadrature points to be used
%        f_scalar - the scalar form of the matrix function f(z) , z scalar
%        f_matrix - A function f_matrix(A,b) computing the action of a matrix function of A 
%        on the vector b 

%Output: fr: Approximation to f(A)b
function fr = rFOM2_v2_invSqrt(p,b,V,H)

if (isempty(p.U))
e1 = zeros(p.m,1);
e1(1)=1;
%Compute approximation
fr = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);
else

term1 = zeros(p.m+p.k,1);

 % Define constant factors appearing in quadrature integration
 Vhat = [p.U V(:,1:p.m)];
 What = [p.C V(:,1:p.m)];
 G = zeros(p.m+p.k,p.m+p.k);
 G(1:p.k,1:p.k) = eye(p.k);
 G(p.k+1:p.m+p.k,p.k+1:p.m+p.k) = H(1:p.m,1:p.m);
 UmC = p.U-p.C;
 e = zeros(p.m,1);
 e(p.m)=1;
 hterm = H(p.m+1,p.m)*V(:,p.m+1)*e';
 VTb = Vhat'*b;
 VTW = Vhat'*What;

 R = @(zx) [zx*UmC  -hterm];
yy = @(zx) (VTW*(zx*speye(p.m+p.k)-G) + Vhat'*R(zx))\VTb;


%compute quadrature nodes and weights
weights = pi/p.num_quad*ones(1,p.num_quad);
t = zeros(1,p.num_quad);
for ii = 1:p.num_quad
t(ii) = cos((2*ii-1)/(2*p.num_quad) * pi);
end
tt = -1*(1-t)./(1+t);

%perform quadrature
for j = 1:p.num_quad
  yterm = yy(tt(j));
  term1 = term1 + weights(j)*(1/(1+t(j)))*yterm;
end

%Compute approximation
fr = (-2/pi)*Vhat*term1;

end