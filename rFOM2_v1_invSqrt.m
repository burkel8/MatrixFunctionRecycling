%rFOM2 - implementation 1

function [deflated_approx] = rFOM2_v1_invSqrt(b,V,H,m,k,U,C,num_quad)
  

term1 = zeros(m,1);
term2 = zeros(k,1);

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


%%Define problem for y
y = @(zx) (zx*speye(m) - H(1:m,1:m) - (zx*VTU - VTC)*( (zx*UTU-UTC)\ ...
    (zx*UTV - UTVp1H) ))\(VTb - (zx*VTU - VTC)*( (zx*UTU - UTC)\UTb));

z1 = @(zx,yx) (zx*UTU - UTC)\(UTb - (zx*UTV - UTVp1H)*yx) ;


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

  yterm = y(tt(j));
  term1 = term1 + (1/(1+t(j)))*weights(1,j)*yterm;
  term2 = term2 + (1/(1+t(j)))*weights(1,j)*z1(tt(j),yterm);
end

%Compute approximation
deflated_approx = const*V(:,1:m)*term1 + const*U*term2;
end