%Function to solve eigenproblem required to compute k  ritz vectors for the matrix A
% using the augmented Krylov subspace km(A,b) + U . 

%Input: m - dimension of Krylov subspace km(A,b)
%       k - number of eigenvectors to compute
%       G - augmented Hessenberg matrix G = [D 0; 0 Hbar];
%       What = [C Vm+1]
%       Vhat = [U Vm+1]

%Output: P - An m+k x k dimensional matrix storing the eigenvectors of the
%small eigenproblem in proposition 8.1 of preprint.
function U = compute_ritz_vectors(p,V,H,order)
    
    P = zeros(p.m+p.k,p.k);

    A = eye(p.m+p.k,p.m+p.k);
    A(1:p.k,p.k+1:end) = p.U'*V(:,1:p.m);
    A(p.k+1:end,1:p.k) = V(:,1:p.m)'*p.U;

    B = [p.U'*p.C  p.U'*(V*H) ;
        V(:,1:p.m)'*p.C H(1:p.m,1:p.m)];

    [harmVecs, harmVals] = eig(A,B);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),order);

    for i=1:p.k
        P(:,i) = harmVecs(:,iperm(i));
    end

    U = [p.U V(:,1:p.m)]*P;
    
end
