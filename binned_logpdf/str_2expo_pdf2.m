    function fvals=str_2expo_pdf2(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        alpha=A(1); a1=A(2); a2=A(3); b=exp(A(4)); lam1=exp(A(5)); lam2=exp(A(6));
        %calculate expo term
        %without xmax
%         stexpo_term=a1*(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))./(gamma_incomplete(b*xmin^alpha,1/alpha));
%         expo_term1=a2*exp(lam1*xmin)*(exp(-lam1*l)-exp(-lam1*u));
%         expo_term2=(1-a1-a2)*exp(lam2*xmin)*(exp(-lam2*l)-exp(-lam2*u));
        %with xmax
        stexpo_term=a1*(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))./(gamma_incomplete(b*xmin^alpha,1/alpha)-gamma_incomplete(b*xmax^alpha,1/alpha));
        expo_term1=a2*(exp(-lam1*l)-exp(-lam1*u))./(exp(-lam1*xmin)-exp(-lam1*xmax));
        expo_term2=(1-a1-a2)*(exp(-lam2*l)-exp(-lam2*u))./(exp(-lam2*xmin)-exp(-lam2*xmax));
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        Pr=log((expo_term2)+(expo_term1)+(stexpo_term));
        fvals=-sum(h.*Pr);
    end