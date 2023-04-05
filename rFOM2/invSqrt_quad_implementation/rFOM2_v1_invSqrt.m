%Function to compute approximation to f(A)b when f is the inverse square root function
% using implementation 1 of r(FOM)^2. The quadrature rule used is outlined
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
function fr = rFOM2_v1_invSqrt(p,b,V,H)

if (isempty(p.U))
e1 = zeros(p.m,1);
e1(1)=1;
%Compute approximation
fr = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);
else

term1 = zeros(p.m,1);
%compute quadrature nodes and weights
weights = pi/p.num_quad*ones(1,p.num_quad);
t = zeros(1,p.num_quad);
for ii = 1:p.num_quad
t(ii) = cos((2*ii-1)/(2*p.num_quad) * pi);
end
tt = -1*(1-t)./(1+t);

term2 = zeros(p.k,1);

%Define constant factors appearing in quadrature integration
VTU = V(:,1:p.m)'*p.U;
VTC = V(:,1:p.m)'*p.C;
UTC = p.U'*p.C;
UTV = p.U'*V(:,1:p.m);
VTb = V(:,1:p.m)'*b;
UTb = p.U'*b;
UTVp1H = p.U'*(V*H);

%function y representing solution to linear system (6.2) in preprint
y = @(zx) (zx*speye(p.m) - H(1:p.m,1:p.m) - (zx*VTU - VTC)*( (zx*eye(p.k)-UTC)\ ...
    (zx*UTV - UTVp1H) ))\(VTb - (zx*VTU - VTC)*( (zx*eye(p.k) - UTC)\UTb));

% function z1 representing the z correction in terms of y
z1 = @(zx,yx) (zx*eye(p.k) - UTC)\(UTb - (zx*UTV - UTVp1H)*yx) ;

%perform quadrature
for j = 1:p.num_quad
  yterm = y(tt(j));
  term1 = term1 + (1/(1+t(j)))*weights(1,j)*yterm;
  term2 = term2 + (1/(1+t(j)))*weights(1,j)*z1(tt(j),yterm);
end

%Compute approximation
fr = (-2/pi)*V(:,1:p.m)*term1 + (-2/pi)*p.U*term2;

end