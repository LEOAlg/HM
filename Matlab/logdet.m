function [ logdetA ] = logDet( A)
    L = chol(A);
logdetA = 2*sum(log(diag(L)));

end

