function  sheet_convert(bra_sheet, gata_sheet, sox_sheet, dist_sheet, plate_id)
    [s,d] = cellfun(@size, bra_sheet);
    mv = max([s,d]);
    bra_m = zeros(length(bra_sheet), mv(2));
    gata_m = zeros(length(bra_sheet), mv(2));
    sox_m = zeros(length(bra_sheet), mv(2));
    dist_m = zeros(length(bra_sheet),mv(2));
    lengths = cellfun(@length, bra_sheet);
    for i = 1:numel(bra_sheet)
        bra_m(i, (1:lengths(i))) = cell2mat(bra_sheet{i});
        sox_m(i, (1:lengths(i))) = cell2mat(sox_sheet{i});
        gata_m(i, (1:lengths(i))) = cell2mat(gata_sheet{i});
        dist_m(i, (1:lengths(i))) = cell2mat(dist_sheet{i});
    end
    writematrix(bra_m,  append('Newbra_sheet_',plate_id,'.csv'));
    writematrix(gata_m, append('Newgata_sheet_',plate_id,'.csv'));
    writematrix(sox_m,  append('Newsox_sheet_',plate_id,'.csv'));
    writematrix(dist_m, append('Newdist_sheet_',plate_id,'.csv'));