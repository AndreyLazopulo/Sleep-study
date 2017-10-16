    function fvals=str_4expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        alpha=A(1); a1=A(2); a2=A(3); a3=A(4); a4=A(5); b=exp(A(6)); lam1=exp(A(7)); lam2=exp(A(8)); lam3=exp(A(9)); lam4=exp(A(10));
        %calculate expo term
        stexpo_term=a1*(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))./(gamma_incomplete(b*xmin^alpha,1/alpha));
        expo_term1=a2*exp(lam1*xmin)*(exp(-lam1*l)-exp(-lam1*u));
        expo_term2=a3*exp(lam2*xmin)*(exp(-lam2*l)-exp(-lam2*u));
        expo_term3=a4*exp(lam3*xmin)*(exp(-lam3*l)-exp(-lam3*u));
        expo_term4=(1-a1-a2-a3-a4)*exp(lam4*xmin)*(exp(-lam4*l)-exp(-lam4*u));
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        Pr=log((expo_term4)+(expo_term3)+(expo_term2)+(expo_term1)+(stexpo_term));
        fvals=-sum(h.*Pr);
    end