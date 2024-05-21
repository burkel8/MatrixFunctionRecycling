function Dkernel = generate_sign_kernel(dim,id)

  fnstr1 = "periodic_L";
  fnstr2 = string(dim);
  fnstr3 = "_b3.55_k0.137n";
  fnstr4 = string(id);
  fnstr5 = ".mat";
  
  fnstr = append(fnstr1,fnstr2,fnstr3,fnstr4,fnstr5);

  load(fnstr);
  Dw = D;
  kappa = 0.137;
  m0 = 1/(2*kappa)-4;
  n = size(Dw,1);
  
  DH = (1/kappa)*(speye(n)-(1/(4+m0))*Dw);
  
  %d = eig(full(DH));
  %plot(real(d),imag(d),'*');
  
  lambda_s = eigs(Dw,1,'smallestreal','Tolerance',1e-10);
  r_s = real(lambda_s);
  m_s = m0 - r_s;
  kappa_c = 1/(2*(4+m_s));

  Dkernel = speye(n) - (4/3)*kappa_c*DH;
  
  %d = eig(full(Dkernel));
  %plot(real(d),imag(d),'*');

end
