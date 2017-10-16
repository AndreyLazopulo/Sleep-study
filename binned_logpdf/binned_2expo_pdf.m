function fvals=binned_2expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        a1=A(1); lam1=exp(A(2)); lam2=exp(A(3));
        %calculate expo term
        expo_term1=a1*exp(lam1*xmin)*(exp(-lam1*l)-exp(-lam1*u));
        expo_term2=(1-a1)*exp(lam2*xmin)*(exp(-lam2*l)-exp(-lam2*u));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        Pr=log((expo_term2)+(expo_term1));
        fvals=-sum(h.*Pr);
    end
    