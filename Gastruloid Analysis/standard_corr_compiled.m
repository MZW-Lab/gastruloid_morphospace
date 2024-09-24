function [g_corr,s_corr,b_corr] = standard_corr_compiled(G_DE_S50, B_DE_S50, S_DE_S50)
    if isempty(G_DE_S50)
        g_corr = 0; 
        s_corr = 0; 
        b_corr = 0;
    else
        g_corr = 1; 
        s_corr = 1; 
        b_corr = 1;
        G_B_corr = corrcoef(G_DE_S50, B_DE_S50);
        G_S_corr = corrcoef(G_DE_S50, S_DE_S50);
        S_B_corr = corrcoef(S_DE_S50, B_DE_S50);
        if G_B_corr == 1
            g_corr = 1; 
        elseif G_S_corr == 1  
            s_corr = 1; 
        elseif S_B_corr == 1       
            b_corr = 1;
        elseif G_B_corr == 0
            g_corr = 0; 
        elseif G_S_corr == 0  
            s_corr = 0; 
        elseif S_B_corr == 0       
            b_corr = 0;
        else
            g_corr = G_B_corr(1, 2); 
            s_corr = G_S_corr(1, 2); 
            b_corr = S_B_corr(1, 2);
        end
    end
    %disp([g_corr, s_corr, b_corr]);
