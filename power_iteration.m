function [lambda_k,b_k] = power_iteration(B)
  %% return largest (in modulus) eigenvalue of inverse of B
  %% B should be equal to (-diff*Laplacian + c - \mu(x))
  
  b_k = rand(size(B,1),1); % random vector
  
  b_k1= zeros(size(B,1),1); err=12; lambda_k1=0; % initialize variables
  
  [L,U,P,Q]=lu(B); % LU decomposition of A
  
  while (err>10^(-8)) 
    lambda_k=lambda_k1; % save previous value
    
    b_k1 = Q*(U\(L\( P*b_k ))); % i.e. b_k1 = A\b_k;
    b_k = b_k1 / norm(b_k1);
    y = Q*(U\(L\( P*b_k ))); % i.e. y=A\b_k;
    lambda_k1 = (y')*(b_k);   
    
    err=abs(lambda_k-lambda_k1);   
  end
end