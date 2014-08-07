function [ return_em, iterator ] = runEM4(Old, W, y_em, i ,epsilon, iteration_limit, y)
    %#codegen
    %% Initilization
    [n,m] = size(y_em);    
    error = Inf;
    likelihood = -Inf;
    I = diag(W);   
    mu = Old.mu;
    C = Old.C;
    sigma = Old.sigma; 
    iterator = 1;  
    pi = 1;
	tau = 1;
	numSamples = sum(W);
    wl = zeros(n,m);
    Cl = zeros(n,n,m);
    while(likelihood < epsilon && iterator < iteration_limit)
    %while(likelihood < epsilon && iterator <= iteration_limit)
        %% E-step for each application  
        Cinv = inv(C);
        temp = inv(eye(n)/sigma + Cinv);  
        
        for l = 1:m            
            if(l~=i)                                                     
                Cl(:,:,l) = temp;
            else              
                Cl(:,:,l) = ( I/sigma + Cinv)^(-1); 
            end  
            wl(:,l)  = Cl(:,:,l) *((y_em(:,l))/sigma + Cinv*mu);                    
        end
        
        %% M-step
        ClSum = (eye(n));
        wlSum = (zeros(n,1));
        wlSumCov = (eye(n));
        normSum = 0;
        for l = 1:m            
            ClSum = ClSum +  Cl(:,:,l);
            wlSum = wlSum + wl(:,l);
            wlSumCov = wlSumCov + (wl(:,l)-mu)*(wl(:,l)-mu)';    
            if(l~=i)                        
            	normSum = normSum + norm(y_em(:,l) - wl(:,l))^2 + trace( Cl(:,:,l)); 
            else
            	normSum = normSum + norm(I*(y_em(:,l) - wl(:,l)))^2 + trace(I* Cl(:,:,l)); 
            end
        end
        mu = ((1/(pi + m))*wlSum);
        C = ((1/(tau + m))*(pi*(mu*mu') + tau*eye(n)+ ClSum + wlSumCov));        
        sigma = ((1/((m-1)*n + numSamples))*normSum);
       % error = norm(Old.mu - mu,'fro')+ norm(Old.C - C,'fro') + Old.sigma - sigma ;    
        Old.mu = mu;        Old.C = C;        Old.sigma = sigma;
        iterator = iterator +1;   
        
        
        %%
        Cinv = inv(C);
        likSum = 0;    
        for l = 1:m            
            if(l~=i)                        
            	likSum = likSum + wl(:,l)'*(eye(n)/sigma + Cinv)*wl(:,l)...
                         -2*(y_em(:,l)'/sigma + mu'*Cinv)*wl(:,l) + y_em(:,l)' *y_em(:,l)/sigma ; 
            
            else
            	likSum = likSum +  wl(:,l)'*(I/sigma + Cinv)*wl(:,l)...
                         -2*(y_em(:,l)'/sigma + mu'*Cinv)*wl(:,l) + y_em(:,l)' *y_em(:,l)/sigma ; 
            end
        end
        likelihood = -(likSum + (sum(W)/m)*log(sigma) -(m+1)*logdet(Cinv) + (m+pi)*mu'*Cinv*mu +trace(Cinv));
%         
        residual = wl(:,i) - y{i};
        rss = sum(residual.*residual); 
        tss = sum((y{i}-mean(y{i})).*(y{i}-mean(y{i})));
        residualsquare = 1-rss/tss;
        accuracy = 1-(1-residualsquare)*(n-1)/(n-2);
       
        fprintf('it = %d, lik = %f, trueAcc = %f\n',iterator,likelihood,accuracy );
%         %fprintf('-------------------------\n');
    	
%        imagesc(Cinv);colorbar;colormap(gray);
%        pause;
    end
    %% Prediction for ith data 
    return_em.w = wl(:,i); 
    return_em.em_t = Old;
    return_em.Cl = Cl;
    return_em.wl = wl; 
end
