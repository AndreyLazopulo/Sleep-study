    function fvals=str_2expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        alpha=A(1); a1=A(2); a2=A(3); b=exp(A(4)); lam1=exp(A(5)); lam2=exp(A(6));
        %calculate expo term
        stexpo_term=log(a1)+log(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))-log(gamma_incomplete(b*xmin^alpha,1/alpha));
        expo_term1=log(a2)+lam1*xmin+log(exp(-lam1*l)-exp(-lam1*u));
        expo_term2=log(1-a1-a2)+lam2*xmin+log(exp(-lam2*l)-exp(-lam2*u));
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        Pr=log(exp(expo_term2)+exp(expo_term1)+exp(stexpo_term));
        fvals=-sum(h.*Pr);
    end