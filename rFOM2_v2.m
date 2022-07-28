%rFOM2 - implementation 2

function [deflated_approx] = rFOM2_v2(b,V,H,m,k,U,C,num_quad, f_scalar)

%Construct contour
s = eigs(H(1:m,1:m),1,'smallestreal');
l = eigs(H(1:m,1:m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(m+k,1);

 % Define constant factors appearing in quadrature integration
 [U,D] = scale_cols_of_U(U,k);
 Vhat = [U V(:,1:m)];
 What = [C V(:,1:m)];
 G = zeros(m+k,m+k);
 G(1:k,1:k) = D;
 G(k+1:m+k,k+1:m+k) = H(1:m,1:m);
 UmC = U-C;
 e = zeros(m,1);
 e(m)=1;
 hterm = H(m+1,m)*V(:,m+1)*e';
 WTb = Vhat'*b;
 WTW = Vhat'*What;

 R = @(zx) [zx*UmC  -hterm];
yy = @(zx) (WTW*(zx*speye(m+k)-G) + Vhat'*R(zx))\WTb;

delta_theta = 2*pi / num_quad;
const = r/num_quad;

%Do quadrature
for j = 1:num_quad

  theta = (j-1)*delta_theta;
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = yy(z);
  term1 = term1 + common_factor*yterm;
  
end

%Compute approximation
deflated_approx = Vhat*const*term1;

end