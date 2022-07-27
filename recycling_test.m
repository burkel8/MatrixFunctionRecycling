%This program runs tests rFOM2 as a recycle method by creating a sequence
%of problems and recycling between them

%Choose parameters for program

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
problem = 'inverse';


%% Parameters of solve
m = 50;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 50;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)
num_quad_points = 1000;   %number of quadrature points (add as many differnt points to this list)
matrix_eps = 0.0;  %parameter to determine how much the matrix changes.
num_systems = 5;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%

e1 = zeros(m,1);
e1(1)=1;

%vectors to store results
err_arnoldi = zeros(1,num_systems);
err_quad_arnoldi = zeros(1,num_systems);
err_rFOM_v1 = zeros(1,num_systems);
err_rFOM_v2 = zeros(1,num_systems);
err_rFOM_v3 = zeros(1,num_systems);
err_recycle_space = zeros(1,num_systems);
eigs_monitor = zeros(1,num_systems);

%store matrix and function in appropriate vectors
[f_scalar, f_matrix] = return_function(problem);
[A,n] = return_matrix(which_matrix,N);

%create vector
b = rand(n,1);
b = b/norm(b);
x = f_matrix(A,b);

%Run Arnoldi on a close matrix to generate first U
Aclose = A + 0.001*sprand(A);
g = rand(n,1);
[Hc,Vc] = arnoldi( Aclose, g , n,m, 1);
[P] = harmonic_ritz(Hc,m,k);
U = Vc(:,1:m)*P;
C = A*U;

%Create sequence of problems and approximate f(A)b for each.
for ix=1:num_systems

    eigs_monitor(ix) = real(eigs(A,1,'smallestreal')); %line not needed

    fprintf("\n Approximating f(A)b # %d ... \n\n",ix);

    %Run Arnoldi
    [H,V] = arnoldi( A, b , n,m, 1);
  
    %Standard Arnoldi Approximation
    arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
    err_arnoldi(ix) = norm(x - arnoldi_approx);

    %Quadrature Arnoldi approximation
    quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad_points,f_scalar);
    err_quad_arnoldi(ix) = norm(x - quad_arnoldi_Approx);

    %rFOM2 f1
    [rFOM_v1_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad_points, f_scalar);
    err_rFOM_v1(ix) = norm(x - rFOM_v1_approx);

    %rFOM2 f2
    [rFOM_v2_approx] = rFOM2_v2(b,V,H,m,k,U,C,num_quad_points, f_scalar);
    err_rFOM_v2(ix) = norm(x - rFOM_v2_approx);

    %rFOM2 f3
    [rFOM_v3_approx] = rFOM2_v3(b,V,H,m,k,U,C,num_quad_points, f_scalar, f_matrix);
    err_rFOM_v3(ix) = norm(x - rFOM_v3_approx);

    fprintf("\n... DONE\n");

     %Construct G
        [U,D] = scale_cols_of_U(U,k);
         Vhat = [U V(:,1:m)];
         What = [C V(:,1:m+1)];
         G = zeros(m+1+k,m+k);
         G(1:k,1:k) = D;
         G(k+1:m+1+k,k+1:m+k) = H;
    
    %Create new problem in sequence and compute its exact solution
    b = rand(n,1);
    b = b/norm(b);
    A = A + matrix_eps*sprand(A);
    x = f_matrix(A,b);

    %Compute new U for next problem in the sequence
    [P] = harm_ritz_aug_krylov(m,k,G,What,Vhat);
    U = Vhat*P;
    [C,R] = qr(A*U,0);
    U = U/R;

         
end

%Plot Results
semilogy(err_arnoldi/norm(x) ,'-s', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
semilogy(err_quad_arnoldi/norm(x) ,'-o', 'LineWidth',1);
hold on;
semilogy(err_rFOM_v1/norm(x),'-v', 'LineWidth',1);
hold on;
semilogy(err_rFOM_v2/norm(x),'-s', 'LineWidth',1,'MarkerSize', 8);
hold on;
semilogy(err_rFOM_v3/norm(x),'-s', 'LineWidth',1);
hold off;

title('sign($\textbf{A}$)\textbf{b} - error vs. problem index','interpreter','latex', 'FontSize', fontsize)
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('$\| f(\textbf{A})\textbf{b} - \textbf{x}_{m} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('Arnoldi','Arnoldi (q)','rFOM$^{2}$ $\tilde{f}_{1}$','rFOM$^{2}$ $\tilde{f}_{2}$','rFOM$^{2}$ $\tilde{f}_{3}$','interpreter','latex');
set(lgd,'FontSize',fontsize);


