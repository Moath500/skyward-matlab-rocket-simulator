function [Coeffs, State] = datcom_parser(mat_name)
linestring = fileread('for006.dat');

%% blocksplit
pattern = '\*+ \<FLIGHT CONDITIONS AND REFERENCE QUANTITIES \*+';
blocks = regexp(linestring,pattern,'split');
blocks = blocks(2:end);

%% get_coeffs_name
pattern =  ' *\<ALPHA\> *([\w -/]*)';
names = {};
for i = 1:length(blocks)
    block = blocks{i};
    token = regexp(block,pattern,'tokens');
   
    % convert cells inside token into strings
    for k = 1:length(token)
        token{k} = char(token{k});
    end
    
    dummy = strjoin(token);
    dummy = split(dummy);
    
    % replaceBadChars
    correct = dummy;
    pattern1 = '[\./-]';
    for j = 1:length(dummy)
        name = dummy{j};
        a = regexprep(name(1:end-1),pattern1,'_');
        b = regexprep(name(end),pattern1,'');
        correct{j} = [a,b];
    end
    
    names{end+1} = correct;
end


%% get_data
pattern1 = ' *\<MACH NO\> *= * (-*[\d]+.[\d]*)';
pattern2 = ' *\<ALTITUDE\> *= * (-*[\d]+.[\d]*)';
pattern3 = ' *\<SIDESLIP\> *= * (-*[\d]+.[\d]*)';
pattern4 = ' *\<MOMENT CENTER\> *= * (-*[\d]+.[\d]*)';

M=[];
A=[];
B=[];
XM=[];

for i = 1:length(blocks)
    block = blocks{i};
    mach = regexp(block,pattern1,'tokens');
    M = [M, str2double(mach{1})];
    sslip = regexp(block,pattern3,'tokens');
    B = [B, str2double(sslip{1})];
    alt = regexp(block,pattern2,'tokens');
    A = [A, str2double(alt{1})];
    mcenter = regexp(block,pattern4,'tokens');
    XM = [XM, str2double(mcenter{1})];
end

%% get_rawdata
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

%% savemat
realM = [M(1)];
realA = [A(1)];
realB = [B(1)];
realXM = [XM(1)];

for i = 2:length(M)
    if not(any(realM == M(i)))
        realM = [realM, M(i)];
    end
    if not(any(realA == A(i)))
        realA = [realA, A(i)];
    end
    if not(any(realB == B(i)))
        realB = [realB, B(i)];
    end
    if not(any(realXM == XM(i)))
        realXM = [realXM, XM(i)];
    end  
end

n_ind = length(realA)*length(realB)*length(realM);
n_rep = length(blocks)/(length(realA)*length(realB)*length(realM));

realNames = {};
for i = 1:n_rep
    realNames = [realNames;names{i}];
end

    for j = 1:length(realNames)
        Coeffs.(realNames{j}) = zeros(length(alpha),length(realM),length(realB),length(realA));
    end

for i = 0:n_ind-1
    index = i*n_rep+1;
    iA = find(realA==A(index));
    iB = find(realB==B(index));
    iM = find(realM==M(index));
    
    for j = 1:length(realNames)
        Coeffs.(realNames{j})(:,iM,iB,iA) = raw_data(:,length(realNames)*i+j);
    end
    
   
end

State.Machs = realM;
State.Alphas = alpha;
State.Betas = realB;
State.Altitudes = realA;

save(mat_name,'State','Coeffs');

end

