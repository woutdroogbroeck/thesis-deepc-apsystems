function [Hp, Hf] = hankel_matrix(col, Tini, Tf)        
        L = Tini+Tf;
    
        H_total = hankel(col(1:L),col(L:end));
        Hp = H_total(1:Tini,:);
        Hf = H_total(Tini+1:end,:);
       
end