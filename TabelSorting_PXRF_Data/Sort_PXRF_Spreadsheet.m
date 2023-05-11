clc; clear; clear all; close all;
fid = fopen("Results March24 2023.csv");
Llines = {};
tline = fgetl(fid);
while ischar(tline)
   Llines{end+1,1} = tline;
   tline = fgetl(fid);
end
fclose(fid);
% Function handle to split each row by comma delimiter 
func = @(input)strsplit(input,',');
% Apply this function to each row in the cell array to get the final data
% As each row is of different lenght we set the UniformOutput to false
raw_data = cellfun(func,Llines,'UniformOutput',false);


var_names = raw_data{1};
N = length(var_names);
M = length(Llines);

header_idx = NaN(M,1);
for i = 1:M
     cl_header = string(raw_data{i,:});
     var_query = string({'File #'});
     if ~isempty(find(strcmp(cl_header,var_query)))
         header_idx(i)=i;
     end
end
header_idx(isnan(header_idx))=[];

% initialize cell array
sorteD = cell(M-length(header_idx),N);

t = 1;
for j = 1:length(header_idx)
    cl_header = raw_data{header_idx(j)};
    k = header_idx(j)+1; % index of the first data after header
    ll = length(raw_data{k});
    lh = length(cl_header);
    
    while(1)
        if ll==lh && isempty(find(ismember(header_idx,k)))
            tmp_data = raw_data{k};
            for s=1:length(cl_header)
                var_query = string(cl_header{s});
                var_list = string(var_names);
                master_header_position=find(strcmp(var_list,var_query));
                if ~isempty(master_header_position)
                    sorteD(t, master_header_position) = tmp_data(s);
                end
            end
            if k~=M % can't exceed length of table
                k = k+1;
                t = t+1;
                tmp_data = raw_data{k};
                ll = length(tmp_data); 
            else
                break
            end
        else
            break
        end
    end
end

writetable(cell2table(sorteD,'VariableNames',string(var_list)),['Sorted Spreadsheet.csv']);

