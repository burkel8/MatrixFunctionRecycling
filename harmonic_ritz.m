%Function to solve harmonic ritz problem
function [P] = harmonic_ritz(H,m,k)

    P = zeros(m,k);
    
    harmRitzMat = speye(m);
    harmRitzMat = harmRitzMat(:,((m-1)+1):m);
    harmRitzMat = H(1:m,:)'\harmRitzMat;
    harmRitzMat = harmRitzMat*(H(m+1:(m+1),(m-1)+1:m)'*H(m+1:(m+1),(m-1)+1:m));
    harmRitzMat = [H(1:m,1:(m-1)) (H(1:m,(m-1)+1:m)+harmRitzMat)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute harmonic Ritz values %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [harmVecs, harmVals] = eig(harmRitzMat);
    harmVals=diag(harmVals);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort values by magnitude %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,iperm] = sort(abs(harmVals));

    idx = 1;
    while(idx <= k)
        P(:,idx) = harmVecs(:,iperm(idx));
        idx = idx + 1;    
    end

end