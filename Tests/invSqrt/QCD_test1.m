% QCD_test1: Tests the quadrature based implementation of r(FOM)2 
% as a recycling method for the special case where the matrix A
% remains fixed and only the right hand sides change. 

% The example used is the sign function of a non-Hermitian QCD matrix of size 3072 x 3072
% Evaluation of the sign function is treated as an evaluation of the inverse square root function.

%The script plots the residual of all approximations as well as the exact error

compute_error = 0; % Binary variable (0 or 1) to decide to also print the exact error
                   % 0 OFF - Do not compute error (Default)
                   % 1 ON  - Compute Error 
order = "descend";

if compute_error == 1
  fprintf('\n Warning: compute_error turned on, code may be slow! \n')
end


addpath(genpath('../'))

%load QCD matrix
load("../../data/smallLQCD.mat");
A = A1;

%Define operator used in the Arnoldi process in this case the operator is A*A
Aop = @(bx) A*(A*bx);


%Define struct p which contains the parameters of the recycling algorithm 
p.n = size(A,1); %problem size
p.m = 100;       % m: Length of Arnoldi cycle
p.k = 30;        % k: Dimension of recycling subspace
p.U = [];        % U: Matrix with columns forming a basis for recycling subspace
p.C = [];        % C: Matrix C given by C = A*U;
p.num_quad = 30; % num_qiuad: Number of quadrature points used in method

% f_scalar: The scalar version of the function f
p.f_scalar = @(zx) 1./sqrt(zx);

%f_matrix: A function which returns the matrix vector product f(A)*b 
p.f_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

num_systems = 30; %Number of differnt f(A)b applications we are testing

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

e1 = zeros(p.m,1);
e1(1)=1;

%vectors to store relative residual for each problem
arnoldi_resid = zeros(1,num_systems);
rfom2_resid = zeros(1,num_systems);

if compute_error == 1
   arnoldi_err = zeros(1,num_systems);
   rfom2_err = zeros(1,num_systems);
end


%Create sequence of problems and approximate f(A)b for each.
for ix=1:num_systems
    fprintf("\n\n\n ### PROBLEM %d ###\n\n\n", ix);
    
    %Create random vector b for f(A)*b application
    b = rand(p.n,1);

    %We actually approximate f(A*A)*(A*b)
    Ab = A*b;

   
    %Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
    [V,H] = arnoldi(Aop, Ab, p);
  
     % Compute Standard Arnoldi Approximation
    fprintf("\n Computing Arnoldi approximation fa...\n");
    fa = norm(Ab)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

    %Compute r(FOM)^{2} approximation
    fprintf("\n Computing r(FOM)^2 approximation fr...\n");
    fr = rFOM2_v1_invSqrt(p,Ab,V,H);


     if compute_error == 1
       exact = p.f_matrix(A*A,Ab);
       arnoldi_err(ix) = norm(exact - fa)/norm(exact);
       rfom2_err(ix) = norm(exact-fr)/norm(exact);
    end

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


    %Estimate Residual of Arnoldi approximation using another call to
    %Arnoldi approximation
    Afa = A*fa;
    [V,H] = arnoldi(Aop,Afa, p);
    faa = norm(Afa)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);
    
    %Estimate residual of r(FOM)^2 approximation using another call to
    %r(FOM)^2
    Afr = A*fr;
    [V,H] = arnoldi(Aop,Afr, p);
    frr = rFOM2_v1_invSqrt(p,Afr,V,H);
   
    %Print Residuals
    arnoldi_resid(ix) = norm(faa - b)/norm(b);
    fprintf("\n Relative Residual (ar) of Arnoldi approximation %2f\n",arnoldi_resid(ix) );

    rfom2_resid(ix) = norm(frr - b)/norm(b);
    fprintf("\n Relative Residual (rr) of r(FOM)^2 approximation %2f\n",rfom2_resid(ix));

    %Pause to give user a chance to read residuals
    pause(5);  

 
end

%Plot Results
if compute_error == 1
semilogy(arnoldi_resid,'-v', 'LineWidth',1);
hold on; 
semilogy(rfom2_resid,'-v', 'LineWidth',1);
hold on;
semilogy(arnoldi_err,'-v', 'LineWidth',1);
hold on; 
semilogy(rfom2_err,'-v', 'LineWidth',1);
hold on;
lgd = legend('Arnoldi residual', 'r(FOM)$^{2}$ residual', 'Arnoldi error', 'r(FOM)$^{2} error$','interpreter','latex');
else
semilogy(arnoldi_resid,'-v', 'LineWidth',1);
hold on; 
semilogy(rfom2_resid,'-v', 'LineWidth',1);
hold off;
lgd = legend('Arnoldi Residual', 'r(FOM)$^{2}$ residual','interpreter','latex');
end
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('', 'FontSize',fontsize);
grid on;

set(lgd,'FontSize',fontsize);


