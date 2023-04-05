% exp_test1: Tests r(FOM)2 as a recycling method for the special case
% where the matrix A remains fixed and only the right hand sides change. 

% The example used is the exponential function of a network matrix of size
% 8297 x 8297

%Previously computed solutions are used to contruct the recycling subspace.
addpath(genpath('../'))

%load network matrix
load("../../data/wiki-Vote.mat");
A = -Problem.A;

%Define operator used in the Arnoldi process in this case the operator is A*A
Aop = @(bx) A*bx;

%Define struct p which contains the parameters of the recycling algorithm 
p.n = size(A,1); %problem size
p.m = 30;       % m: Length of Arnoldi cycle
p.k = 10;        % k: Dimension of recycling subspace
p.U = [];        % U: Matrix with columns forming a basis for recycling subspace
p.C = [];        % C: Matrix C given by C = A*U;
p.num_quad = 500; % num_qiuad: Number of quadrature points used in method

% f_scalar: The scalar version of the function f
p.f_scalar = @(zx) exp(zx);

%f_matrix: A function which returns the matrix vector product f(A)*b 
p.f_matrix = @(Ax,bx) expm(full(Ax))*bx;

num_systems = 5; %Number of differnt f(A)b applications we are testing

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

e1 = zeros(p.m,1);
e1(1)=1;

%vectors to store relative residual for each problem
arnoldi_err = zeros(1,num_systems);
rfom2_err = zeros(1,num_systems);


%Create sequence of problems and approximate f(A)b for each.
for ix=1:num_systems
    fprintf("\n\n\n ### PROBLEM %d ###\n\n\n", ix);
    
    %Create random vector b for f(A)*b application
    b = rand(p.n,1);

    exact = p.f_matrix(A,b);

    %Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
    [V,H] = arnoldi(Aop, b, p);
  
     % Compute Standard Arnoldi Approximation
    fprintf("\n Computing Arnoldi approximation fa...\n");
    fa = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

    %Compute r(FOM)^{2} approximation
    fprintf("\n Computing r(FOM)^2 approximation fr...\n");
    fr = rFOM2_v3(p,A,b,V,H);

     %update U and c
    fprintf("\n Updating U and C ... \n")

    %Add newly computed approximate solution to recycling subspace.
    U(:,ix) = fr;
    p.k = size(U,2);
    [p.U,~] = qr(U,0);
    p.C = Aop(p.U);

    arnoldi_err(ix) = norm(exact - fa)/norm(exact);
    fprintf("\n Relative error (er) of arnoldi approximation %2f\n",arnoldi_err(ix));
    rfom2_err(ix) = norm(exact - fr)/norm(exact);
   fprintf("\n Relative error (er) of r(FOM)^2 approximation %2f\n",rfom2_err(ix));
    
end

%Plot Results
semilogy(arnoldi_err,'-v', 'LineWidth',1);
hold on; 
semilogy(rfom2_err,'-v', 'LineWidth',1);
hold off;
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('Relative Error', 'FontSize',fontsize);
grid on;
lgd = legend('Arnoldi Error', 'r(FOM)$^{2}$ Error','interpreter','latex');
set(lgd,'FontSize',fontsize);


