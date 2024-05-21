%This function computes the standard Arnoldi approximation using for f(A)b
%where f is the inverse square root function. A brief outline of the
%quadrature rule is given in section 7 of the preprint and references
%therin.

%Inputs: 
%        V,H - basis of Krylov subspace and hessenberg matrix H constructed
%        from Arnoldi process to build basis for Km(A,b)
%        m - dimension of Krylov subspace
%        num_quad - number of quadrature points to be used

%Output: fa: Approximation to f(A)b
function fa = quad_arnoldi_invSqrt(p,b,V,H)

e1 = zeros(p.m,1);
e1(1) = 1;
term1 = zeros(p.m,1);

%Function y representing the FOM approximation to the shifted linear system
%(sig*I - A)x(sig) = b
y = @(zx) (zx*speye(p.m)-H(1:p.m,1:p.m))\e1;


%Compute quadrature nodes and weights
weights = pi/p.num_quad*ones(1,p.num_quad);
t = zeros(1,p.num_quad);
for ii = 1:p.num_quad
t(ii) = cos((2*ii-1)/(2*p.num_quad) * pi);
end
tt = -1*(1-t)./(1+t);

%perform quadrature
for j = 1:p.num_quad
   term1 = term1 + weights(j)*(1/(1+t(j)))*y(tt(j));
end

%Compute approximation
fa = norm(b)*(-2/pi)*V(:,1:p.m)*term1;
end