function fvals=binned_5expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        a1=A(1); a2=A(2); a3=A(3); a4=A(4); lam1=exp(A(5)); lam2=exp(A(6)); lam3=exp(A(7)); lam4=exp(A(8));  lam5=exp(A(9));
        %calculate expo term
        expo_term1=a1*exp(lam1*xmin)*(exp(-lam1*l)-exp(-lam1*u));
        expo_term2=a2*exp(lam2*xmin)*(exp(-lam2*l)-exp(-lam2*u));
        expo_term3=a3*exp(lam3*xmin)*(exp(-lam3*l)-exp(-lam3*u));
        expo_term4=a4*exp(lam4*xmin)*(exp(-lam4*l)-exp(-lam4*u));
        expo_term5=(1-a1-a2-a3-a4)*exp(lam5*xmin)*(exp(-lam5*l)-exp(-lam5*u));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        expo_term5 = reshape(expo_term5, 1, numel(expo_term5));
        Pr=log((expo_term5)+(expo_term4)+(expo_term3)+(expo_term2)+(expo_term1));
        fvals=-sum(h.*Pr);
    end