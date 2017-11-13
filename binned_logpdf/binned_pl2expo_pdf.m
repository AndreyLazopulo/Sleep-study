function fvals=binned_pl2expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        alpha=A(1); a1=A(2); lam1=exp(A(3)); lam2=exp(A(4));
        %calculate expo term
        %without xmax
%         expo_term1=a1*(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))./(gamma_incomplete(lam1*xmin,1-alpha));
%         expo_term2=(1-a1)*(gamma_incomplete(lam2*l,1-alpha)-gamma_incomplete(lam2*u,1-alpha))./(gamma_incomplete(lam2*xmin,1-alpha));
        %with xmax
        expo_term1=a1*(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))./(gamma_incomplete(lam1*xmin,1-alpha)-gamma_incomplete(lam1*xmax,1-alpha));
        expo_term2=(1-a1)*(gamma_incomplete(lam2*l,1-alpha)-gamma_incomplete(lam2*u,1-alpha))./(gamma_incomplete(lam2*xmin,1-alpha)-gamma_incomplete(lam2*xmax,1-alpha));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        Pr=log((expo_term1)+(expo_term2));
        fvals=-sum(h.*Pr);
    end