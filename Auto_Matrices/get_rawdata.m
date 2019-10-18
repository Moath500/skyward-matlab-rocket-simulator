function [alpha,raw_data] = get_rawdata(blocks)

pattern = '^[-\d](?=[\d\.])';
pattern2 = '\n\t?';

raw_data = {};

for i = 1:length(blocks)
    block = blocks{i};
    lines = regexp(block,pattern2,'split');
    new_1 = {};
    
    for j = 1:length(lines)
        line = lines{j};
        line = strip(line);
        
        if not(isempty(regexp(line,pattern)))
            new_1{end+1} = str2double(split(line));
        end
        
    end
    
    alpha = [];
    new_1_1 = {};
    
    for j = 1:length(new_1)
        row = new_1{j};
        alpha = [alpha, row(1)];
        new_1_1{end+1} = row(2:end);
    end
    new_1 = new_1_1;
    l = length(new_1{1});
    
    
    for j = 1:length(new_1_1)
        row = new_1_1{j};
        if length(row)~=l
            new_1{j-length(new_1_1)/2} = [ new_1{j-length(new_1_1)/2}; new_1_1{j}];
            new_1{j} = {};
        else
            index = j;
        end
    end
    
    new_1 = transpose(cell2mat(new_1(1:index)));
    
    raw_data{end+1} = new_1;
            
end

if index ~= length(new_1_1)
    alpha = alpha(1:end/2);
end

raw_data = cell2mat(raw_data);

end