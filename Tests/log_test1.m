% log_test1: Tests the quadrature based implementation of r(FOM)2 
% as a recycling method for the special case where the matrix A
% remains fixed and only the right hand sides change. 

% The example used is the log function of a Poisson matrix

order = "descend";

addpath(genpath('../'))

set(0,...
 'defaultaxeslinewidth',1,...
'defaultaxesfontsize',18,...
'defaultlinelinewidth',2,...
'defaultpatchlinewidth',2,...
'defaultlinemarkersize',8,...
'defaulttextinterpreter','latex');

A = gallery("poisson",40);
n = size(A,1);

Aop = @(bx) A*bx;

%Define struct p which contains the parameters of the recycling algorithm 
p.n = size(A,1); %problem size
p.m = 50;       % m: Length of Arnoldi cycle
p.k = 20;        % k: Dimension of recycling subspace
p.U = [];        % U: Matrix with columns forming a basis for recycling subspace
p.C = [];        % C: Matrix C given by C = A*U;
p.num_quad = 30000; % num_qiuad: Number of quadrature points used in method

% f_scalar: The scalar version of the function f
p.f_scalar = @(zx) log(zx);

%f_matrix: A function which returns the matrix vector product f(A)*b 
p.f_matrix = @(Ax,bx) logm(full(Ax))*bx;

num_systems = 20; %Number of differnt f(A)b applications we are testing

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

e1 = zeros(p.m,1);
e1(1)=1;

arnoldi_err = zeros(1,num_systems);
rfom2_err = zeros(1,num_systems);

%Create sequence of problems and approximate f(A)b for each.
for ix=1:num_systems
    fprintf("\n\n\n ### PROBLEM %d ###\n\n\n", ix);
    
    %Create random vector b for f(A)*b application
    b = rand(p.n,1);

    %Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
    [V,H] = arnoldi(Aop, b, p);
  
     % Compute Standard Arnoldi Approximation
    fprintf("\n Computing Arnoldi approximation fa...\n");
    fa = norm(b)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

    %Compute r(FOM)^{2} approximation
    fprintf("\n Computing r(FOM) approximation fr...\n");
    fr = rFOM2_v1(p,b,V,H);

    exact = p.f_matrix(A,b);
    arnoldi_err(ix) = norm(exact - fa)/norm(exact);
    rfom2_err(ix) = norm(exact-fr)/norm(exact);
   
    %update U and c
    fprintf("\n Updating U and C ... \n")

    if (isempty(p.U))
    [P,~] = eigs(H(1:p.m,1:p.m),p.k,'smallestabs');
    U = V(:,1:p.m)*P;
    else
    U = compute_harmonic_ritz_vectors(p,V,H,order);
    end

    [p.U,~] = qr(U,0);
    p.C = Aop(p.U);

end


semilogy(arnoldi_err,'-v');
hold on; 
semilogy(rfom2_err,'-v');
hold off;
lgd = legend('FOM', 'rFOM','interpreter','latex');
xlabel('Problem index');
ylabel("Relative error");
grid on;



