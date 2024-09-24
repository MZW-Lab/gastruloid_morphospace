function mah_result = mahalanobis(gata_sheet, bra_sheet, sox_sheet)

p_ctrl_I = find(strcmp(drug_name, 'BMP4_1'));
n_ctrl_I = find(strcmp(drug_name, 'mTeSR1_1'));
p_ctrl_GMI = cell2mat(gata_sheet(p_ctrl_I, :));
n_ctrl_GMI = cell2mat(gata_sheet(n_ctrl_I, :));

p_ctrl_BMI = cell2mat(bra_sheet(p_ctrl_I, :));
n_ctrl_BMI = cell2mat(bra_sheet(n_ctrl_I, :));

p_ctrl_SMI = cell2mat(sox_sheet(p_ctrl_I, :));
n_ctrl_SMI = cell2mat(sox_sheet(n_ctrl_I, :));

ctrl_length = length(p_ctrl_SMI);

for i = 1:length(drug_name)
    gata_y = cell2mat(gata_sheet(i, :));
    bra_y = cell2mat(bra_sheet(i, :));
    sox_y = cell2mat(sox_sheet(i, :));
    p_ctrl_Gadj = p_ctrl_GMI;
    n_ctrl_Gadj = n_ctrl_GMI;
    p_ctrl_Badj = p_ctrl_BMI;
    n_ctrl_Badj = n_ctrl_BMI;
    p_ctrl_Sadj = p_ctrl_SMI;
    n_ctrl_Sadj = n_ctrl_SMI;
    y_cell_len = length(sox_y);
    if ctrl_length > y_cell_len
        p_ctrl_Gadj = p_ctrl_GMI(1:y_cell_len);
        n_ctrl_Gadj = n_ctrl_GMI(1:y_cell_len);
        p_ctrl_Badj = p_ctrl_BMI(1:y_cell_len);
        n_ctrl_Badj = n_ctrl_BMI(1:y_cell_len);
        p_ctrl_Sadj = p_ctrl_SMI(1:y_cell_len);
        n_ctrl_Sadj = n_ctrl_SMI(1:y_cell_len);
    elseif ctrl_length < y_cell_len
        gata_y = gata_y(1:ctrl_length);
        bra_y = bra_y(1:ctrl_length);
        sox_y  = sox_y(1:ctrl_length);
    end
    pc_d_mahal_gata = mahal(gata_y,p_ctrl_Gadj); 
    nc_d_mahal_gata = mahal(gata_y,n_ctrl_Gadj);
    pc_d_mahal_bra = mahal(bra_y,p_ctrl_Badj);
    nc_d_mahal_bra = mahal(bra_y,n_ctrl_Badj);
    pc_d_mahal_sox = mahal(sox_y,p_ctrl_Sadj);
    nc_d_mahal_sox = mahal(sox_y,n_ctrl_Sadj);
end
pc_d_mahal_gata = mahal(n_ctrl_GMI,p_ctrl_GMI); 
pc_d_mahal_bra = mahal(n_ctrl_BMI,p_ctrl_BMI);
pc_d_mahal_sox = mahal(n_ctrl_SMI,p_ctrl_SMI);

scatter3(p_ctrl_GMI,p_ctrl_SMI ,p_ctrl_BMI);
hold on
scatter3(n_ctrl_GMI,n_ctrl_SMI ,n_ctrl_BMI);
hold offm
scatter3(pc_d_mahal_gata,pc_d_mahal_sox,pc_d_mahal_bra, 'o');

pc_d_mahal_gata_f = mahal(p_ctrl_GMI,n_ctrl_GMI); 
pc_d_mahal_bra_f = mahal(p_ctrl_BMI,n_ctrl_BMI);
pc_d_mahal_sox_f = mahal(p_ctrl_SMI,n_ctrl_SMI);
