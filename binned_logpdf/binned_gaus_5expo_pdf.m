function fvals=binned_gaus_5expo_pdf(x,bin_bound,A)
        [h,edges]=histcounts(x,bin_bound);
        h = reshape(h, 1, numel(h));
        edges = reshape(edges, 1, numel(edges));
        l=edges(1:end-1);
        u=edges(2:end);
        xmin=min(x);
        xmax=max(x);
        a1=A(1); a2=A(2); a3=A(3); a4=A(4); a5=A(5); mu=A(6); sig=A(7); 
        lam2=exp(A(8)); lam3=exp(A(9)); lam4=exp(A(10)); lam5=exp(A(11)); lam6=exp(A(12));
        %calculate expo term
        %without xmax
%         expo_term1=a1*exp(lam1*xmin)*(exp(-lam1*l)-exp(-lam1*u));
%         expo_term2=a2*exp(lam2*xmin)*(exp(-lam2*l)-exp(-lam2*u));
%         expo_term3=a3*exp(lam3*xmin)*(exp(-lam3*l)-exp(-lam3*u));
%         expo_term4=(1-a1-a2-a3)*exp(lam4*xmin)*(exp(-lam4*l)-exp(-lam4*u));
        %with xmax
        gaussian_term=a1*(erf((mu-u)/(sqrt(2)*sig))-erf((mu-l)/(sqrt(2)*sig)))./(erf((mu-xmax)/(sqrt(2)*sig))-erf((mu-xmin)/(sqrt(2)*sig))); 
        expo_term2=a2*(exp(-lam2*l)-exp(-lam2*u))./(exp(-lam2*xmin)-exp(-lam2*xmax));
        expo_term3=a3*(exp(-lam3*l)-exp(-lam3*u))./(exp(-lam3*xmin)-exp(-lam3*xmax));
        expo_term4=a4*(exp(-lam4*l)-exp(-lam4*u))./(exp(-lam4*xmin)-exp(-lam4*xmax));
        expo_term5=a5*(exp(-lam5*l)-exp(-lam5*u))./(exp(-lam5*xmin)-exp(-lam5*xmax));
        expo_term6=(1-a1-a2-a3-a4-a5)*(exp(-lam6*l)-exp(-lam6*u))./(exp(-lam6*xmin)-exp(-lam6*xmax));
        gaussian_term = reshape(gaussian_term, 1, numel(gaussian_term));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        expo_term5 = reshape(expo_term5, 1, numel(expo_term5));
        expo_term6 = reshape(expo_term6, 1, numel(expo_term6));
        Pr=log((expo_term6)+(expo_term5)+(expo_term4)+(expo_term3)+(expo_term2)+(gaussian_term));
        fvals=-sum(h.*Pr);
    end