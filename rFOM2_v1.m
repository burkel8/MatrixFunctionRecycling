%rFOM2 - implementation 1

function [deflated_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad, f_scalar)

%Construct contour
s = eigs(H(1:m,1:m),1,'smallestreal');
l = eigs(H(1:m,1:m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(m,1);
term2 = zeros(k,1);
term3 = zeros(k,1);

%Define constant factors appearing in quadrature integration
%Denote Vp = V(:,1:m+1) and V = V(:,1:m) and a "T" for transpose
VTU = V(:,1:m)'*U;
VTC = V(:,1:m)'*C;
UTU = U'*U;
UTC = U'*C;
UTV = U'*V(:,1:m);
VTb = V(:,1:m)'*b;
UTb = U'*b;
UTVp1H = U'*V*H;


%%Define functions (both definitions of y below are equivelant)

%y = @(zx) (...
%           zx*speye(m) - H(1:m,1:m) ...
%           - zx*(zx*VTU - VTC)*((zx*UTU - UTC)\UTV)...
%           +(zx*VTU - VTC)*((zx*UTU - UTC)\UTVp1H))\ ...
%           (VTb - (zx*VTU - VTC)*((zx*UTU - UTC)\UTb) ) ;

y = @(zx) (zx*speye(m) - H(1:m,1:m) - (zx*VTU - VTC)*( (zx*UTU-UTC)\ ...
    (zx*UTV - UTVp1H) ))\(VTb - (zx*VTU - VTC)*( (zx*UTU - UTC)\UTb));

z1 = @(zx,yx) (zx*UTU - UTC)\(UTb - (zx*UTV - UTVp1H)*yx) ;

%Do quadrature
delta_theta = 2*pi / num_quad;
const = r/num_quad;
for j = 1:num_quad

  theta = (j-1)*delta_theta;
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = y(z);
  term1 = term1 + common_factor*yterm;
  term2 = term2 + common_factor*z1(z,yterm);
end

%Compute approximation
deflated_approx = const*V(:,1:m)*term1 + const*U*term2;
end