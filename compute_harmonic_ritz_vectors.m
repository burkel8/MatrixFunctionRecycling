%Function to solve eigenproblem required to compute k harmonic ritz vectors for the matrix A
% using the augmented Krylov subspace km(A,b) + U . 

%Input: m - dimension of Krylov subspace km(A,b)
%       k - number of eigenvectors to compute
%       G - augmented Hessenberg matrix G = [D 0; 0 Hbar];
%       What = [C Vm+1]
%       Vhat = [U Vm+1]

%Output: P - An m+k x k dimensional matrix storing the eigenvectors of the
%small eigenproblem in proposition 8.1 of preprint.
function U = compute_harmonic_ritz_vectors(p,V,H,order)
    
    P = zeros(p.m+p.k,p.k);

    %Construct G
     Vhat = [p.U V(:,1:p.m)];
    What = [p.C V(:,1:p.m+1)];
    G = zeros(p.m+1+p.k,p.m+p.k);
    G(1:p.k,1:p.k) = eye(p.k);
    G(p.k+1:p.m+1+p.k,p.k+1:p.m+p.k) = H;



    B = G'*(What'*What)*G;
    A = G'*(What'*Vhat);



    % A = [p.C'*p.U  p.C'*V(:,1:p.m) ;
    %     H'*(V'*p.U)   H(1:p.m,1:p.m)];
    % 
    % B = [p.C'*p.C  p.C'*(V*H) ;
    %     H'*(V'*p.C) H'*H];

    [harmVecs, harmVals] = eig(A,B);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),order);

    for i=1:p.k
        P(:,i) = harmVecs(:,iperm(i));
    end

    U = Vhat*P;
    
end

