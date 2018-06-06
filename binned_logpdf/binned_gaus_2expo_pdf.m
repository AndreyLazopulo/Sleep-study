function fvals=binned_gaus_2expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        a1=A(1); a2=A(2); mu=A(3); sig=A(4); lam2=exp(A(5)); lam3=exp(A(6));
        %calculate expo term
        %without xmax
%         expo_term1=a1*exp(lam1*xmin)*(exp(-lam1*l)-exp(-lam1*u));
%         expo_term2=a2*exp(lam2*xmin)*(exp(-lam2*l)-exp(-lam2*u));
%         expo_term3=a3*exp(lam3*xmin)*(exp(-lam3*l)-exp(-lam3*u));
%         expo_term4=(1-a1-a2-a3)*exp(lam4*xmin)*(exp(-lam4*l)-exp(-lam4*u));
        %with xmax
        gaussian_term=a1*(erf((mu-u)/(sqrt(2)*sig))-erf((mu-l)/(sqrt(2)*sig)))./(erf((mu-xmax)/(sqrt(2)*sig))-erf((mu-xmin)/(sqrt(2)*sig))); 
        expo_term2=a2*(exp(-lam2*l)-exp(-lam2*u))./(exp(-lam2*xmin)-exp(-lam2*xmax));
        expo_term3=(1-a1-a2)*(exp(-lam3*l)-exp(-lam3*u))./(exp(-lam3*xmin)-exp(-lam3*xmax));
        gaussian_term = reshape(gaussian_term, 1, numel(gaussian_term));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        Pr=log((expo_term3)+(expo_term2)+(gaussian_term));
        fvals=-sum(h.*Pr);
    end