which_matrix = "smallLQCD";
problem = "inverse";
m = 40;
k = 20;
num_quad_points = 5000;
num_systems = 10;
e1 = zeros(m,1);
e1(1)=1;

err_arnoldi = zeros(1,num_systems);
err_quad_arnoldi = zeros(1,num_systems);
err_rFOM = zeros(1,num_systems);
err_rGMRES = zeros(1,num_systems);
err_recycle_space = zeros(num_systems);
% loading the system matrix
if which_matrix=="large_sch"
    load mat_schwinger128x128b3phasenum11000.mat;
    n = size(S,1);
    shift = 0.0;
elseif which_matrix=="small_sch"
    load Schw_16.mat
    n = size(S,1);
    S(:,n/2+1:n) = -S(:,n/2+1:n); % ---> do this for the 16^{2} Schwinger
    shift = 0.0;
elseif which_matrix == "smallLQCD"
      load('smallLQCD_A1.mat')
      S=A1;
      shift = 0;
      [n,~] = size(S); 
elseif which_matrix == "poisson"

N = 50;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
n = N*N;
s = eigs(A,1,'smallestabs');
S = A;
else
   fprintf('ERR: no matrix chosen!');
end

%S = S + shift*speye(n);

f1_scalar = @(zx) 1.0/zx;
f1_matrix = @(Ax,bx) Ax\bx;

f2_scalar = @(zx) 1.0/sqrt(zx);
f2_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

if problem =="inverse"
    f_scalar = @(zx) f1_scalar(zx);
    f_matrix = @(Ax,bx) f1_matrix(Ax,bx);
elseif problem == "sign"
    f_scalar = @(zx) f2_scalar(zx);
    f_matrix = @(Ax,bx) f2_matrix(Ax,bx);
else
    error("ERROR : unknown function chosen!\n");
end


b = rand(n,1);
b = b/norm(b);
x = f_matrix(S,b);
x0 = zeros(n,1);

%% Run Arnoldi on a close matrix
Sclose = S + 0.001*sprand(S);
g = rand(n,1);
[Hc,Vc] = arnoldi( Sclose, g , n,m, 1);
[P] = harmonic_ritz(Hc,m,k);
U = Vc(:,1:m)*P;
C = S*U;


theta = (U'*U)\(U'*S*U);
eigs_res_norm = norm(U*theta - S*U)/norm(U);
fprintf("avg_res_norm = %f\n", eigs_res_norm);


%% computing f(A)b

for ix=1:num_systems

    fprintf("\nSOLVING FOR SYSTEM # %d ... \n\n",ix);
    defl_subsp_tol = 1.0e-3;
    [H,V] = arnoldi( S, b , n,m, 1);
  
     arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
     err_arnoldi(ix) = norm(x - arnoldi_approx);


    quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad_points,f_scalar);
    err_quad_arnoldi(ix) = norm(x - quad_arnoldi_Approx);

    [rFOM_approx] = rFOM2(b,V,H,m,k,U,C,num_quad_points, f_scalar);
    err_rFOM(ix) = norm(x - rFOM_approx);

    [rGMRES_approx] = rGMRES2(b,V,H,m,k,U,C,num_quad_points, f_scalar);
    err_rGMRES(ix) = norm(x - rGMRES_approx);

    fprintf("\n... DONE\n");

       [U,D] = scale_cols_of_U(U,k);
        Vhat = [U V(:,1:m)];
        What = [C V(:,1:m+1)];
        G = zeros(m+1+k,m+k);
        G(1:k,1:k) = D;
      
        G(k+1:m+1+k,k+1:m+k) = H;
        [P] = harm_ritz_aug_krylov(m,k,G,What,Vhat);
        Y = Vhat*P;
        [Q,R]=qr(G*P,0);
        C = What*Q;
        U = Y/R;


 theta = (U'*U)\(U'*S*U);
eigs_res_norm = norm(U*theta - S*U)/norm(U);
fprintf("avg_res_norm = %f\n", eigs_res_norm);
err_recycle_space(ix) = eigs_res_norm;


    b = rand(n,1);
    b = b/norm(b);
    %S = S + 0.0001*sprand(S);
    x = f_matrix(S,b);
    
end

xx=1:1:num_systems;
semilogy(xx,err_arnoldi/norm(x) ,'-x');
hold on;
semilogy(xx,err_quad_arnoldi/norm(x) ,'-s');
hold on;
semilogy(xx,err_rFOM/norm(x),'-o');
hold on;
semilogy(xx,err_rGMRES/norm(x),'-o');
hold off;

title('sign($\textbf{A}$)\textbf{b} - approximation accuracy vs number of quadrature nodes','interpreter','latex')
xlabel('System number','interpreter','latex');
ylabel('$\| f(A)b - x_{m} \|_{fro}$','interpreter','latex');
grid on;
legend('Arnoldi','Arnoldi via quadrature','rFOM$^{2}$','rGMRES$^{2}$','interpreter','latex');
xticks(xx);

