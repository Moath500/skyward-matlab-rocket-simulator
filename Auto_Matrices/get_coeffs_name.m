function cnames = get_coeffs_name(blocks)
pattern =  ' *\<ALPHA\> *([\w -/]*)';
cnames = {};
for i = 1:length(blocks)
    block = blocks{i};
    token = regexp(block,pattern,'tokens');
    % convert cells inside token into strings
    for k = 1:length(token)
        token{k} = char(token{k});
    end
    dummy = strjoin(token);
    dummy = split(dummy);
    dummy = replaceBadChars(dummy);
    cnames{end+1} = dummy;    
end
