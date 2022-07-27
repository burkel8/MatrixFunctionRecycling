%This file runs the same code experiment as in recycling_test.m for the
%sign function of a lattice QCD matrix corresponding to an 8^{4} lattice.
% The exact solutions were precomputed and are loaded accordingly 

load("data/exact_sign_LQCD_8to4_0.mat");
load("data/LQCD_8to4_0.mat");
load("data/rhs_sign_LQCD_8to4_0.mat");
[n,~] = size(A);
A1 = A;

%User parameters
m = 50;
k = 20;
num_quad_points = 1000;
num_systems = 5;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

e1 = zeros(m,1);
e1(1)=1;

%vectors to store results
err_arnoldi = zeros(1,num_systems);
err_quad_arnoldi = zeros(1,num_systems);
err_rFOM_v1 = zeros(1,num_systems);
err_rFOM_v2 = zeros(1,num_systems);
err_rFOM_v3 = zeros(1,num_systems);

f2_scalar = @(zx) 1.0/sqrt(zx);
f2_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

f_scalar = @(zx) f2_scalar(zx);
f_matrix = @(Ax,bx) f2_matrix(Ax,bx);

% Run Arnoldi on a close matrix to construct first U
Aclose = A + 0.001*sprand(A);
g = rand(n,1);
[Hc,Vc] = arnoldi( Aclose, g , n,m, 1);
[P] = harmonic_ritz(Hc,m,k);
U = Vc(:,1:m)*P;
C = A*U;

%computing f(A)b for each system.
for ix=1:num_systems

    fprintf("\n Approximating f(A)b # %d ... \n\n",ix);
    
    %Run Arnoldi 
    [H,V] = arnoldi( A, b , n,m, 1);
  
    %Arnoldi approximation
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
 
    %Construct variables for harmritz problem
    [U,D] = scale_cols_of_U(U,k);
    Vhat = [U V(:,1:m)];
    What = [C V(:,1:m+1)];
    G = zeros(m+1+k,m+k);
    G(1:k,1:k) = D;
    G(k+1:m+1+k,k+1:m+k) = H;

    %load next problem
     if ix < num_systems
        load(['data/exact_sign_LQCD_8to4_', num2str(ix), '.mat']);
        load(['data/LQCD_8to4_', num2str(ix), '.mat']);
        load(['data/rhs_sign_LQCD_8to4_', num2str(ix), '.mat']);
     end
 
     %construct U for next system.
     [P] = harm_ritz_aug_krylov(m,k,G,What,Vhat);
     U = Vhat*P;
     [C,R] = qr(A*U,0);
     U = U/R;
end

%plot results
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




