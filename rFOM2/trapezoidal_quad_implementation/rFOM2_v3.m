%Function to compute approximation to f(A)b using implementation 3 of
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

function fr = rFOM2_v3(p,A,b,V,H)

if (isempty(p.U))
e1 = zeros(p.m,1);
e1(1)=1;
%Compute approximation
fr = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

else

%Construct circular contour with radius r and centre circle_centre such
%that the spectrum of A is contained within the contour.
%s = eigs(H(1:p.m,1:p.m),1,'smallestreal');
%l = eigs(H(1:p.m,1:p.m),1,'largestreal');

s = eigs(A,1,'smallestreal');
l = eigs(A,1,'largestreal');

sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(p.m+p.k,1);

 % Define constant factors appearing in quadrature integration
 Vhat = [p.U V(:,1:p.m)];
 What = [p.C V(:,1:p.m)];
 G = zeros(p.m+p.k,p.m+p.k);
 G(1:p.k,1:p.k) = eye(p.k);
 G(p.k+1:p.m+p.k,p.k+1:p.m+p.k) = H(1:p.m,1:p.m);
 UmC = p.U-p.C;

em = zeros(p.m,1); em(p.m) = 1.0;
VTb = Vhat'*b;
VTW = Vhat'*What;
VTWinvVTb = (VTW)\VTb;

hterm = -H(p.m+1,p.m)*V(:,p.m+1)*em';
I = speye(p.m+p.k);

R = @(zx) [zx*UmC,hterm];
Gz = @(zx) VTW*(zx*speye(p.k+p.m)-G) ;
B = @(zx) Vhat'*R(zx);

%Function yy representing the quadrature integral I in eq (6.8) of preprint.
yy = @(zx) Gz(zx)\((I + B(zx)*(Gz(zx)\I))\(B(zx)*(Gz(zx)\VTb)));

delta_theta = 2*pi / p.num_quad;

%perform quadrature
for j = 1:p.num_quad

  theta = (j-1)*delta_theta;

  %move to new point in circle.
  z = r*exp(1i*theta) + circle_centre;
  common_factor = p.f_scalar(z)*exp(1i*theta);
  yterm = yy(z);
  term1 = term1 + common_factor*yterm;
end

%Compute Approximation. First term is closed form term in eq (6.8) in
%preprint. Second term comes from the quadrature.
fr = Vhat*p.f_matrix(G,VTWinvVTb) - (r/p.num_quad)*Vhat*term1;
end