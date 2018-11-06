function [mylogpdfvals]=my4ExpoLogPdf(x,p,l,k,xmax,xmin)
        Ap=[p,l];
        lam=(Ap(2*k-1)); a=Ap(1:k-1); 
        
        %%%last term
        constant = 1-exp(-lam*(xmax-xmin));
        terms = (1-sum(a))*lam.*exp(-lam*(x-xmin))./constant;
        %%% rest terms
        for i=1:k-1
            a=Ap(i); lam=(Ap(k+i-1));
            constant=1-exp(-lam*(xmax-xmin));
            newterm=a*lam.*exp(-lam.*(x-xmin))./constant;
            terms=terms+newterm;
        end        
        mylogpdfvals=terms;
    end