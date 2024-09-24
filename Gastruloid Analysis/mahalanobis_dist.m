function mal_info = mahalanobis_dist(gata_sheet, bra_sheet, sox_sheet, drug_name, )

n_ctrl_I = find(strcmp(drug_name, 'mtser1'));
p_ctrl_I = find(strcmp(drug_name, 'bmp4'));

n_ctrl_gX = cell2mat(gata_sheet{n_ctrl_I,:});  %#ok<FNDSB>
p_ctrl_gX = cell2mat(gata_sheet{p_ctrl_I,:});  %#ok<FNDSB>
n_ctrl_bX = cell2mat(bra_sheet{n_ctrl_I,:});
p_ctrl_bX = cell2mat(bra_sheet{p_ctrl_I,:});
n_ctrl_sX = cell2mat(sox_sheet{n_ctrl_I,:});
p_ctrl_sX = cell2mat(sox_sheet{p_ctrl_I,:});
n_ctrl_len = length(n_ctrl_gX);
p_ctrl_len = length(p_ctrl_gX);
for i = 1:length(drug_name)
    gata_y = cell2mat(gata_sheet{i,:});
    bra_y = cell2mat(bra_sheet{i,:});
    sox_y = cell2mat(sox_sheet{i,:});
    y_len = length(gata_y);
    if y_len > n_ctrl_len
       gata_y_n = gata_y(1:n_ctrl_len);
       bra_y_n = bra_y(1:n_ctrl_len);
       sox_y_n = sox_y(1:n_ctrl_len);
    else
       gata_y_n = gata_y;
       bra_y_n = bra_y;
       sox_y_n = sox_y;
    end
    if y_len > p_ctrl_len
       gata_y_p = gata_y(1:p_ctrl_len);
       bra_y_p = bra_y(1:p_ctrl_len);
       sox_y_p = sox_y(1:p_ctrl_len);
    else
       gata_y_p = gata_y;
       bra_y_p = bra_y;
       sox_y_p = sox_y;
    end
   gata_p_mahal = mahal(gata_y_p', p_ctrl_gX');
   gata_n_mahal = mahal(gata_y_n', n_ctrl_gX');
   bra_p_mahal = mahal(bra_y_p', p_ctrl_bX');
   bra_n_mahal = mahal(bra_y_n', n_ctrl_bX');
   sox_p_mahal = mahal(sox_y_p', p_ctrl_sX');
   sox_n_mahal = mahal(sox_y_n', n_ctrl_sX');
end