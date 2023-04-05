%Function to compute approximation to f(A)b using implementation 1 of
%r(FOM)^2. The quadrature rule used is a trapezoidal rule.

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

function fr = rFOM2_v1(p,b,V,H)


if (isempty(p.U))
e1 = zeros(p.m,1);
e1(1)=1;
%Compute approximation
fr = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);
else

%Construct circular contour with radius r and centre circle_centre such
%that the spectrum of A is contained within the contour.
s = eigs(H(1:p.m,1:p.m),1,'smallestreal');
l = eigs(H(1:p.m,1:p.m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(p.m,1);
term2 = zeros(p.k,1);

%Define constant factors appearing in quadrature integration
VTU = V(:,1:p.m)'*p.U;
VTC = V(:,1:p.m)'*p.C;
UTC = p.U'*p.C;
UTV = p.U'*V(:,1:p.m);
VTb = V(:,1:p.m)'*b;
UTb = p.U'*b;
UTVp1H = p.U'*V*H;

%function y representing solution to linear system (6.2) in preprint
y = @(zx) (zx*speye(p.m) - H(1:p.m,1:p.m) - (zx*VTU - VTC)*( (zx*eye(p.k)-UTC)\ ...
    (zx*UTV - UTVp1H) ))\(VTb - (zx*VTU - VTC)*( (zx*eye(p.k) - UTC)\UTb));

% function z1 representing the z correction in terms of y
z1 = @(zx,yx) (zx*eye(p.k) - UTC)\(UTb - (zx*UTV - UTVp1H)*yx) ;

delta_theta = 2*pi / p.num_quad;

%Perform quadrature
for j = 1:p.num_quad

  theta = (j-1)*delta_theta;

  %move to new quadrature point
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = p.f_scalar(z)*exp(1i*theta);
  yterm = y(z);
  term1 = term1 + common_factor*yterm;
  term2 = term2 + common_factor*z1(z,yterm);
end

%Compute approximation
fr = (r/p.num_quad)*V(:,1:p.m)*term1 + (r/p.num_quad)*p.U*term2;

end