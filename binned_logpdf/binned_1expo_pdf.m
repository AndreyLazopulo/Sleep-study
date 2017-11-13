function fvals=binned_1expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        a1=1; lam1=exp(A(1));
        %calculate expo term
%         expo_term1=log(a1)+lam1*xmin+log(exp(-lam1*l)-exp(-lam1*u));
        expo_term1=log(a1)+log(exp(-lam1*l)-exp(-lam1*u))-log(exp(-lam1*xmin)-exp(-lam1*xmax));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        Pr=expo_term1;
        fvals=-sum(h.*Pr);
    end
    