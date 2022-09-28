%Function which returns the matrix to be used in each experiment.
%User can add a matrix here if they wish.
function [A,n] = return_matrix(which_matrix,N,mass)

if which_matrix == "smallLQCD"
      load('smallLQCD_A1.mat');
      A=A1;
      [n,~] = size(A); 
elseif which_matrix == "hermetian_QCD"
      load("hermetian_QCD.mat");
      A = Problem.A;
      n = size(A,2);
     
      eyel = eye(256);
      eye3 = eye(3);
      gamma5 = [0 0 -1 0; 0 0 0 -1; -1 0 0 0; 0 -1 0 0];
      k1 = kron(eyel,gamma5);
      Gamma5 = kron(k1,eye3);
      Q = Gamma5*A;
      Q2 = Q'*Q;
      A = Q2;
      Q = Q - mass*speye(n);
      A = Q;
elseif which_matrix == "poisson"
 
A = (N+1)^2*gallery('poisson',N);
n = N*N;

elseif which_matrix == "chemical_potential"
D2 = (N+1)^2*gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = D2 + 0*D1;
S = A;
[n,~]=size(S);
else    
    fprintf('ERR: no matrix chosen!');
end

end


