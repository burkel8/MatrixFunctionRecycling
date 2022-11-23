function [standard_approx] = quad_arnoldi_invSqrt(V,H,m,num_quad)

e1 = zeros(m,1);
e1(1) = 1;
term1 = zeros(m,1);


y = @(zx) (zx*speye(m)-H(1:m,1:m))\e1;

const = -2/pi;

 weights = pi/num_quad*ones(1,num_quad);
 %The notes then are the points themselves.
%In this case we take them between zero and one
t = zeros(1,num_quad);
for ii = 1:num_quad
t(ii) = cos((2*ii-1)/(2*num_quad) * pi);
end
 tt = -1*(1-t)./(1+t);



for j = 1:num_quad
 
  term1 = term1 + weights(j)*(1/(1+t(j)))*y(tt(j));


end

standard_approx = const*V(:,1:m)*term1;

end