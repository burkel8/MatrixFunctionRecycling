function [deflated_approx] = rGMRES2(b,V,H,m,k,U,C,num_quad,f_scalar)

 
s = eigs(H(1:m,1:m),1,'smallestabs');
l = eigs(H(1:m,1:m),1,'largestabs');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     



term1 = zeros(m,1);
term2 = zeros(k,1);
term3 = zeros(k,1);
%% Define constant factors appearing in quadrature integration

%Denote Vp = V(:,1:m+1) and V = V(:,1:m) and a "T" for transpose
VpTV = V(:,1:m+1)'*V(:,1:m);
HTH = H'*H;
VpTU = V(:,1:m+1)'*U;
VpTC = V(:,1:m+1)'*C;
CTU = C'*U;
CTC = C'*C;
CTVpH = C'*V(:,1:m+1)*H;
CTb = C'*b;
HTVp1Tb = H'*V(:,1:m+1)'*b;
CTV = C'*V(:,1:m);

y = @(zx) (...
         zx*H'*VpTV - ...
        HTH -... 
       zx*H'*(zx*VpTU - VpTC)*((zx*CTU - CTC)\CTV) + ...
        H'*(zx*VpTU - VpTC)*((zx*CTU - CTC)\ CTVpH ))\ ...
       (HTVp1Tb - H'*(zx*VpTU - VpTC)*((zx*CTU - CTC)\ CTb ));

z1 = @(zx) ((zx*CTU - CTC)\CTb);

z2 = @(zx) (zx*CTU - CTC)\(zx*CTV - CTVpH);

delta_theta = 2*pi / num_quad;
const = (1i*r*delta_theta)/(2*pi*1i);
[n,~] = size(V);
for j = 1:num_quad
  theta = (j-1)*delta_theta;
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = y(z);
  term1 = term1 + common_factor*yterm;
  term2 = term2 + common_factor*z1(z);
  term3 = term3 + common_factor*z2(z)*yterm;


%   proj = eye(n) - (z*U-C)*((z*C'*U - C'*C)\C');
%   rank(V(:,1:m))
%   rank(proj'*(V*H))

end

deflated_approx = const*V(:,1:m)*term1 + const*U*term2 - const*U*term3;

end