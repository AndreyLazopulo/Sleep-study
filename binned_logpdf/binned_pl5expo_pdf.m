function fvals=binned_pl5expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        alpha=A(1); a1=A(2); a2=A(3); a3=A(4); a4=A(5); lam1=exp(A(6)); lam2=exp(A(7)); lam3=exp(A(8)); lam4=exp(A(9)); lam5=exp(A(10));
        %calculate expo term
        expo_term1=a1*(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))./(gamma_incomplete(lam1*xmin,1-alpha));
        expo_term2=a2*(gamma_incomplete(lam2*l,1-alpha)-gamma_incomplete(lam2*u,1-alpha))./(gamma_incomplete(lam2*xmin,1-alpha));
        expo_term3=a3*(gamma_incomplete(lam3*l,1-alpha)-gamma_incomplete(lam3*u,1-alpha))./(gamma_incomplete(lam3*xmin,1-alpha));
        expo_term4=a4*(gamma_incomplete(lam4*l,1-alpha)-gamma_incomplete(lam4*u,1-alpha))./(gamma_incomplete(lam4*xmin,1-alpha));
        expo_term5=(1-a1-a2-a3-a4)*(gamma_incomplete(lam5*l,1-alpha)-gamma_incomplete(lam5*u,1-alpha))./(gamma_incomplete(lam5*xmin,1-alpha));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        expo_term5 = reshape(expo_term5, 1, numel(expo_term5));
        Pr=log((expo_term1)+(expo_term2)+(expo_term3)+(expo_term4)+(expo_term5));
        fvals=-sum(h.*Pr);
    end