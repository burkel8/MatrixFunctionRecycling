% test_changin_mat: Tests the quadrature based implementation of r(FOM)
% as a recycling method for the special case where the matrix A
% is slowly changing.

% The example used is the sign function of a non-Hermitian QCD matrix of size 3072 x 3072
% Evaluation of the sign function is treated as an evaluation of the inverse square root function.

%The script plots the residual of all approximations as well as the exact err1or
addpath(genpath('../'))

order = "descend";

set(0,...
 'defaultaxeslinewidth',1,...
'defaultaxesfontsize',18,...
'defaultlinelinewidth',2,...
'defaultpatchlinewidth',2,...
'defaultlinemarkersize',8,...
'defaulttextinterpreter','latex');

%Define struct p which contains the parameters of the recycling algorithm
n = 3072;       %problem size
p.n=n;
p.m = 50;       % m: Length of Arnoldi cycle
p.k = 20;        % k: Dimension of recycling subspace
p.U = [];        % U: Matrix with columns forming a basis for recycling subspace
p.C = [];        % C: Matrix C given by C = A*U;
p.num_quad = 30; % num_qiuad: Number of quadrature points used in method
epsilon = 0;

% f_scalar: The scalar version of the function f
p.f_scalar = @(zx) 1./sqrt(zx);

%f_matrix: A function which returns the matrix vector product f(A)*b
p.f_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

num_systems = 20; %Number of differnt f(A)b applications we are testing

e1 = zeros(p.m,1);
e1(1)=1;

arnoldi_err = zeros(1,num_systems);
rfom2_err2 = zeros(1,num_systems);

%Create sequence of problems and approximate f(A)b for each.
for ix=1:num_systems
    
    fprintf("\n\n\n ### PROBLEM %d ###\n\n\n", ix);

    load(['data/4to4/periodic_L4_b3.55_k0.137n0_' num2str(ix) '.mat'])
    A = D;

    %% The real deal
    Aop = @(bx) A'*(A*bx);

    if ~isempty(p.U)
       p.C = Aop(p.U);
    end

    sqA = sqrtm(full(A'*A));
    eigs(A'*A,1,'smallestreal')

    %Create random vector b for f(A)*b application
    b = rand(p.n,1);

    Ab = b;
    %Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
    [V,H] = arnoldi(Aop, Ab, p);

    % Compute Standard Arnoldi Approximation
    fprintf("\n Computing Arnoldi approximation...\n");
    fa = norm(Ab)*V(:,1:p.m)*p.f_matrix(H(1:p.m,1:p.m),e1);

    %Compute r(FOM)^{2} approximation
    fprintf("\n Computing r(FOM)^2 v2 approximation...\n");
    frr = rFOM2_v2_invSqrt(p,Ab,V,H);

    exact = sqA\Ab;
    arnoldi_err(ix) = norm(exact - fa)/norm(exact);
    rfom2_err2(ix) = norm(exact-frr)/norm(exact);


    %update U and c
    fprintf("\n Updating U and C ... \n")

    if (isempty(p.U))
    
        [P,~] = eigs(H(1:p.m,1:p.m),p.k,'smallestabs');
        U = V(:,1:p.m)*P;

    else

        U = compute_ritz_vectors(p,V,H,order);
    end

     [p.U,~] = qr(U,0);

end


semilogy(arnoldi_err,'--');
hold on;
semilogy(rfom2_err2,':v');
hold on;
lgd = legend( 'FOM', 'rFOM','interpreter','latex');
grid on;
xlabel("Problem index");
ylabel("Relative error");
set(lgd);
