function fvals=binned_pl1expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        alpha=A(1); a1=1; lam1=exp(A(2));
        %calculate expo term
        expo_term1=log(a1)+log(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))-log(gamma_incomplete(lam1*xmin,1-alpha));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        Pr=expo_term1;
        fvals=-sum(h.*Pr);
    end