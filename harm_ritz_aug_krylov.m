%Function to solve augmented harmonic ritz problem.
function [P] = harm_ritz_aug_krylov(m, k, G, What, Vhat)
    
    P = zeros(m+k,k);

    %Remember for our case columns of C and V are no longer orthonormal so we 
    % do not have B = G'*G;. Instead we have the following

    B = G'*(What'*What)*G;
    A = G'*(What'*Vhat);

    [harmVecs, harmVals] = eig(A,B);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),'descend');

    for i=1:k
        P(:,i) = harmVecs(:,iperm(i));
    end
end

