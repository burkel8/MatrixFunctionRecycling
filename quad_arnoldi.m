function [standard_approx] = quad_arnoldi(b,V,H,m,num_quad, f_scalar)

e1 = zeros(m,1); e1(1)=1;

s = eigs(H(1:m,1:m),1,'smallestreal');
l = eigs(H(1:m,1:m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(m,1);
y = @(zx) (zx*speye(m)-H(1:m,1:m))\e1;

delta_theta = 2*pi / num_quad;
const = (1i*r*delta_theta)/(2*pi*1i);

for j = 1:num_quad
  theta = (j-1)*delta_theta;
  z = r*exp(1i*theta) + circle_centre;
  common_factor = f_scalar(z)*exp(1i*theta);
  term1 = term1 + common_factor*y(z);
end

standard_approx = norm(b)*const*V(:,1:m)*term1;

end