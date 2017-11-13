function fvals=binned_pl6expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        alpha=A(1); a1=A(2); a2=A(3); a3=A(4); a4=A(5); a5=A(6); 
        lam1=exp(A(7)); lam2=exp(A(8)); lam3=exp(A(9)); lam4=exp(A(10)); lam5=exp(A(11)); lam6=exp(A(12));
        %calculate expo term
        %without xmax
%         expo_term1=a1*(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))./(gamma_incomplete(lam1*xmin,1-alpha));
%         expo_term2=a2*(gamma_incomplete(lam2*l,1-alpha)-gamma_incomplete(lam2*u,1-alpha))./(gamma_incomplete(lam2*xmin,1-alpha));
%         expo_term3=a3*(gamma_incomplete(lam3*l,1-alpha)-gamma_incomplete(lam3*u,1-alpha))./(gamma_incomplete(lam3*xmin,1-alpha));
%         expo_term4=a4*(gamma_incomplete(lam4*l,1-alpha)-gamma_incomplete(lam4*u,1-alpha))./(gamma_incomplete(lam4*xmin,1-alpha));
%         expo_term5=a5*(gamma_incomplete(lam5*l,1-alpha)-gamma_incomplete(lam5*u,1-alpha))./(gamma_incomplete(lam5*xmin,1-alpha));
%         expo_term6=(1-a1-a2-a3-a4-a5)*(gamma_incomplete(lam6*l,1-alpha)-gamma_incomplete(lam6*u,1-alpha))./(gamma_incomplete(lam6*xmin,1-alpha));
        %with xmax
        expo_term1=a1*(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))./(gamma_incomplete(lam1*xmin,1-alpha)-gamma_incomplete(lam1*xmax,1-alpha));
        expo_term2=a2*(gamma_incomplete(lam2*l,1-alpha)-gamma_incomplete(lam2*u,1-alpha))./(gamma_incomplete(lam2*xmin,1-alpha)-gamma_incomplete(lam2*xmax,1-alpha));
        expo_term3=a3*(gamma_incomplete(lam3*l,1-alpha)-gamma_incomplete(lam3*u,1-alpha))./(gamma_incomplete(lam3*xmin,1-alpha)-gamma_incomplete(lam3*xmax,1-alpha));
        expo_term4=a4*(gamma_incomplete(lam4*l,1-alpha)-gamma_incomplete(lam4*u,1-alpha))./(gamma_incomplete(lam4*xmin,1-alpha)-gamma_incomplete(lam4*xmax,1-alpha));
        expo_term5=a5*(gamma_incomplete(lam5*l,1-alpha)-gamma_incomplete(lam5*u,1-alpha))./(gamma_incomplete(lam5*xmin,1-alpha)-gamma_incomplete(lam5*xmax,1-alpha));
        expo_term6=(1-a1-a2-a3-a4-a5)*(gamma_incomplete(lam6*l,1-alpha)-gamma_incomplete(lam6*u,1-alpha))./(gamma_incomplete(lam6*xmin,1-alpha)-gamma_incomplete(lam6*xmax,1-alpha));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        expo_term5 = reshape(expo_term5, 1, numel(expo_term5));
        expo_term6 = reshape(expo_term6, 1, numel(expo_term6));
        Pr=log((expo_term1)+(expo_term2)+(expo_term3)+(expo_term4)+(expo_term5)+(expo_term6));
        fvals=-sum(h.*Pr);
    end