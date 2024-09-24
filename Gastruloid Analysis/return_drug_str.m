%% Function Description
%pulls the names of each respective drug file and stores them as file naming agents

%% Return Drug String
function [samplenames, drug_name, drug_str, drug_str_pdf, filenames_mod] = return_drug_str( filename )

samplenames = [];
filenames_mod = filename(arrayfun(@(x) x.name(1), filename) ~= ('.')); %Remove '.','..' from array call
toRemove = strcmp({filenames_mod.name}, 'Thumbs.db');
filenames_mod(toRemove) = [];
%splits the names of the files for access
for ii = 1:numel(filenames_mod)
    tmp = strsplit(filenames_mod(ii).name, {'_','.tif'});
    drug_name{ii} = strrep(tmp{1},' ', '_');
    drug_str{ii} = append(drug_name{ii},'.fig');
    drug_str_pdf{ii} = append(drug_name{ii},'.pdf');
    samplenames{ii} = tmp{end-1};
end