    function fvals=str_1expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        alpha=A(1); a1=A(2); b=exp(A(3)); lam=exp(A(4));
        %calculate expo term
        %without xmax
%         stexpo_term=a1*(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))./(gamma_incomplete(b*xmin^alpha,1/alpha));
%         expo_term=(1-a1)*exp(lam*xmin)*(exp(-lam*l)-exp(-lam*u));
        %with xmax
        stexpo_term=a1*(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))./(gamma_incomplete(b*xmin^alpha,1/alpha)-gamma_incomplete(b*xmax^alpha,1/alpha));
        expo_term=(1-a1)*(exp(-lam*l)-exp(-lam*u))./(exp(-lam*xmin)-exp(-lam*xmax));
        
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        expo_term = reshape(expo_term, 1, numel(expo_term));
        Pr=log((expo_term)+(stexpo_term));
        fvals=-sum(h.*Pr);
    end