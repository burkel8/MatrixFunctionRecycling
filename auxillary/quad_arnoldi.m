%This function computes the standard Arnoldi approximation using
%quadrature. The quadrature rule used here is the trapezoidal rule.

%Inputs: b - vector b for which we want to approximate f(A)b
%        V,H - basis of Krylov subspace and hessenberg matrix H constructed
%        from Arnoldi process to build basis for Km(A,b)
%        m - dimension of Krylov subspace
%        num_quad - number of quadrature points to be used
%        f_scalar - the scalar form of the matrix function f(z) , z scalar

%Output: fa: Approximation to f(A)b
function fa = quad_arnoldi(p)

e1 = zeros(p.m,1); e1(1)=1;
%set up a circular contour with radius r and centre circle_centre.
% Ensure spectrum of A lies within contour.
s = eigs(p.H(1:p.m,1:p.m),1,'smallestreal');
l = eigs(p.H(1:p.m,1:p.m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2;
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(p.m,1);

%Function y representing shifted FOM approximation to shifted linear system.
y = @(zx) (zx*speye(p.m)-p.H(1:p.m,1:p.m))\e1;

delta_theta = 2*pi / p.num_quad; %distance between quadrature points
const = (1i*r*delta_theta)/(2*pi*1i);

%Perform quadrature
for j = 1:p.num_quad
  theta = (j-1)*delta_theta;
  z = r*exp(1i*theta) + circle_centre;
  common_factor = p.f_scalar(z)*exp(1i*theta);
  term1 = term1 + common_factor*y(z);
end

%compute final approximation
fa = norm(p.b)*const*p.V(:,1:p.m)*term1;
end