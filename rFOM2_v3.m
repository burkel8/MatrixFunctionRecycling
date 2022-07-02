function [deflated_approx] = rFOM2_v3(b,V,H,m,k,U,C,num_quad, f_scalar, f_matrix)

s = eigs(H(1:m,1:m),1,'smallestreal');
l = eigs(H(1:m,1:m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(m+k,1);

 %% Define constant factors appearing in quadrature integration
 [U,D] = scale_cols_of_U(U,k);
 Vhat = [U V(:,1:m)];
 What = [C V(:,1:m)];
 G = zeros(m+k,m+k);
 G(1:k,1:k) = D;
 G(k+1:m+k,k+1:m+k) = H(1:m,1:m);
 UmC = U-C;

em = zeros(m,1); em(m) = 1.0;
WTb = What'*b;
WTW = What'*What;
WTWinvWTb = (WTW)\WTb;

hterm = -H(m+1,m)*V(:,m+1)*em';
I = speye(m+k);

R = @(zx) [zx*UmC,hterm];
Gz = @(zx) WTW*(zx*speye(k+m)-G) ;
B = @(zx) What'*R(zx);

%We found it more numerically stable to keep WTb in the function yy
%instead of having it in the final update once.
yy = @(zx) Gz(zx)\((I + B(zx)*(Gz(zx)\I))\(B(zx)*(Gz(zx)\WTb)));

delta_theta = 2*pi / num_quad;
const = r/num_quad;

for j = 1:num_quad

  theta = (j-1)*delta_theta;
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = yy(z);
  term1 = term1 + common_factor*yterm;
  
end


deflated_approx = Vhat*f_matrix(G,WTWinvWTb) - const*Vhat*term1;

end