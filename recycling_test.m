%% Step 1: Choose parameters for program

%%First choose the matrix possible options are
%-- A small lattice QCD matrix of size 3072x3072 ("smallLQCD")
%-- A poisson matrix of size N*N x N*N (user specifies N) ("poisson")
%-- A chemical potential matrix of size N*N x N*N (user specifies N) ("chemical_potantial")
which_matrix = "smallLQCD";   

%%Choose the function . Available options are
% -- inverse function ("inverse")
% -- Sign function ("sign")
% -- log function ("log")
% -- square root function ("sqrt")
problem = 'sign';


%% Parameters of solve
m = 50;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 50;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)
num_quad_points = 10000;   %number of quadrature points (add as many differnt points to this list)
matrix_eps = 0.001;  %parameter to determine how much the matrix changes.
num_systems = 10;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%

e1 = zeros(m,1);
e1(1)=1;

err_arnoldi = zeros(1,num_systems);
err_quad_arnoldi = zeros(1,num_systems);
err_rFOM_v1 = zeros(1,num_systems);
err_rFOM_v2 = zeros(1,num_systems);
err_rFOM_v3 = zeros(1,num_systems);
err_recycle_space = zeros(1,num_systems);
eigs_monitor = zeros(1,num_systems);

[f_scalar, f_matrix] = return_function(problem);
[S1,n] = return_matrix(which_matrix,N);
S = S1;

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

    eigs_monitor(ix) = real(eigs(S,1,'smallestreal'));

    fprintf("\nSOLVING FOR SYSTEM # %d ... \n\n",ix);
    defl_subsp_tol = 1.0e-3;
    [H,V] = arnoldi( S, b , n,m, 1);
  
     arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
     err_arnoldi(ix) = norm(x - arnoldi_approx);


    quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad_points,f_scalar);
    err_quad_arnoldi(ix) = norm(x - quad_arnoldi_Approx);

    [rFOM_v1_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad_points, f_scalar);
    err_rFOM_v1(ix) = norm(x - rFOM_v1_approx);

    [rFOM_v2_approx] = rFOM2_v2(b,V,H,m,k,U,C,num_quad_points, f_scalar);
    err_rFOM_v2(ix) = norm(x - rFOM_v2_approx);

    [rFOM_v3_approx] = rFOM2_v3(b,V,H,m,k,U,C,num_quad_points, f_scalar, f_matrix);
    err_rFOM_v3(ix) = norm(x - rFOM_v3_approx);

   fprintf("\n... DONE\n");

        [U,D] = scale_cols_of_U(U,k);
         Vhat = [U V(:,1:m)];
         What = [C V(:,1:m+1)];
         G = zeros(m+1+k,m+k);
         G(1:k,1:k) = D;
         G(k+1:m+1+k,k+1:m+k) = H;
        
        [P] = harm_ritz_aug_krylov(m,k,G,What,Vhat);
        U = Vhat*P;




    b = rand(n,1);
    b = b/norm(b);
    S = S + matrix_eps*sprand(S);
    x = f_matrix(S,b);

    %Maake U and C compaitable with new system
    [Q,R]=qr(S*U,0);
     C = Q;
     U = U/R;

      theta = (U'*U)\(U'*S1*U);
eigs_res_norm = norm(U*theta - S1*U)/norm(U);
fprintf("avg_res_norm = %f\n", eigs_res_norm);
err_recycle_space(ix) = eigs_res_norm;
    
end

xx=1:1:num_systems;

semilogy(xx,err_arnoldi/norm(x) ,'-s', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
semilogy(xx,err_quad_arnoldi/norm(x) ,'-o', 'LineWidth',1);
hold on;
semilogy(xx,err_rFOM_v1/norm(x),'-v', 'LineWidth',1);
hold on;
semilogy(xx,err_rFOM_v2/norm(x),'-s', 'LineWidth',1,'MarkerSize', 8);
hold on;
semilogy(xx,err_rFOM_v3/norm(x),'-s', 'LineWidth',1);
hold off;

title('sign($\textbf{A}$)\textbf{b} - error vs. problem index','interpreter','latex', 'FontSize', fontsize)
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('$\| f(\textbf{A})\textbf{b} - \textbf{x}_{m} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('Arnoldi','Arnoldi (q)','rFOM$^{2}$ $\tilde{f}_{1}$','rFOM$^{2}$ $\tilde{f}_{2}$','rFOM$^{2}$ $\tilde{f}_{3}$','interpreter','latex');
set(lgd,'FontSize',fontsize);
xticks(xx);


