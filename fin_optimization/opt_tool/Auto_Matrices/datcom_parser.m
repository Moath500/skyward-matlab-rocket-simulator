function [Coeffs, State] = datcom_parser(mat_name)
linestring = fileread('for006.dat');

%% blocksplit
pattern = '\*+ \<FLIGHT CONDITIONS AND REFERENCE QUANTITIES \*+';
blocks = regexp(linestring,pattern,'split');
blocks = blocks(2:end);

%% get_coeffs_name
pattern =  ' *\<ALPHA\> *([\w -/]*)';
names = cell(26,1);
index = 1;

for i = 1:4
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
    
    names(index:index+length(correct)-1) = correct;
    index = index+length(correct);
end


%% get_data
pattern1 = ' *\<MACH NO\> *= * (-*[\d]+.[\d]*)';
pattern2 = ' *\<ALTITUDE\> *= * (-*[\d]+.[\d]*)';
pattern3 = ' *\<SIDESLIP\> *= * (-*[\d]+.[\d]*)';
pattern4 = ' *\<MOMENT CENTER\> *= * (-*[\d]+.[\d]*)';

M=zeros(1,length(blocks)/4);
A=zeros(1,length(blocks)/4);
B=zeros(1,length(blocks)/4);
XM=zeros(1,length(blocks)/4);

for i = 1:length(blocks)/4
    block = blocks{(i-1)*4+1};
    mach = regexp(block,pattern1,'tokens');
    M(i) = str2double(mach{1});
    sslip = regexp(block,pattern3,'tokens');
    B(i) = str2double(sslip{1});
    alt = regexp(block,pattern2,'tokens');
    A(i) = str2double(alt{1});
    mcenter = regexp(block,pattern4,'tokens');
    XM(i) = str2double(mcenter{1});
end

%% get_alpha
pattern = '^[-\d](?=[\d\.])';
pattern2 = '\n\t?';

block = blocks{2}; 
lines = regexp(block,pattern2,'split');
index = 0;
new_1 = cell(200,1);

for j = 1:length(lines)
    line = lines{j};
    line = strip(line);
    
    if regexp(line,pattern,'once')
        index = index + 1;
        new_1{index} = str2double(split(line));
    end
    
end

alpha = zeros(1,index);

for j = 1:index
    row = new_1{j};
    alpha(j) = row(1);
end
%% get_rawdata
raw_data = cell(1,length(blocks));

for i = 1:length(blocks)
    block = blocks{i};
    lines = regexp(block,pattern2,'split');
    index = 0;
    new_1 = cell(200,1);
    
    for j = 1:length(lines)
        line = lines{j};
        line = strip(line);
        
        if regexp(line,pattern,'once')
            index = index + 1;
            new_1{index} = str2double(split(line));
        end
        
    end
    
    new_1_1 = cell(1,index);
    
    for j = 1:index
        row = new_1{j};
        new_1_1{j} = row(2:end);
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
    
    raw_data{i} = new_1;
            
end


raw_data = cell2mat(raw_data);

%% savemat
realM = [M(1), zeros(1,200)];
realA = [A(1), zeros(1,200)];
realB = [B(1), zeros(1,200)];

iM = 1;
iA = 1;
iB = 1;

for i = 2:length(M)
    if not(any(realM == M(i)))
        iM = iM + 1;
        realM(iM) = M(i);
    end
    if not(any(realA == A(i)))
        iA = iA + 1;
        realA(iA) = A(i);
    end
    if not(any(realB == B(i)))
        iB = iB + 1;
        realB(iB) = B(i);
    end
end

realM = realM(1:iM);
realA = realA(1:iA);
realB = realB(1:iB);


    for j = 1:length(names)
        Coeffs.(names{j}) = zeros(length(alpha),iM,iB,iA);
    end

for i = 1:length(blocks)/4
    index = i;
    iA = find(realA==A(index));
    iB = find(realB==B(index));
    iM = find(realM==M(index));
    
    for j = 1:length(names)
        Coeffs.(names{j})(:,iM,iB,iA) = raw_data(:,length(names)*(i-1)+j);
    end
    
   
end

State.Machs = realM;
State.Alphas = alpha;
State.Betas = realB;
State.Altitudes = realA;

save(mat_name,'State','Coeffs');

end

