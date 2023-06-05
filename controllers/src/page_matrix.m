function [Pp, Pf] = page_matrix(col, Tini, Tf, sigma)        
        L = Tini+Tf;
        T = length(col);
        
        total_col = T/L;
        P_total = zeros(L, total_col);
%         for j = 1:T-L+1
%             P_total(:,j) = col((j-1)*L+1 : j*L);
%         end
        k = 1;
        
        for j = 1:total_col
            for i = 1:L
                P_total(i,j) = col(k);
                k = k +1;
            end
        end
        if ~ isempty(sigma)
            [U,S,V] = svd(P_total);   
            S(S<sigma)=0;
            P_total = U*S*V';
        end
        Pp = P_total(1:Tini,:);
        Pf = P_total(Tini+1:end,:);
       
end