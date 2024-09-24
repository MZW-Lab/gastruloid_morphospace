i1 = []; i2 = []; i3 = [];

for i = 1:25
    inds = find(cell_dist_bins25==i);
    fprintf('%i\n', i);
    v1 = sum([sub_struct.(col_name).C2(inds).MeanIntensity]);
    v2 = sum([sub_struct.(col_name).C3(inds).MeanIntensity]);
    v3 = sum([sub_struct.(col_name).C4(inds).MeanIntensity]);
    i1(i) = v1;
    i2(i) = v2;
    i3(i) = v3;
    
end 
av1 = i1./cell_dist_binc25;
av2 = i2./cell_dist_binc25;
av3 = i3./cell_dist_binc25;