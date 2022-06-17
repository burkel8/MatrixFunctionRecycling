function [P] = harmonic_ritz(H,m,k)

    P = zeros(m,k);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build Matrix for Eigenvalue problem %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select k smallest eigenvectors                                 %
    % Optionally store k+1 vectors to capture complex conjugate pair %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx = 1;
    while(idx <= k)
        P(:,idx) = harmVecs(:,iperm(idx));
        idx = idx + 1;    
    end

end