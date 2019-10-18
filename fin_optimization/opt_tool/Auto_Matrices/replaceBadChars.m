function correct = replaceBadChars(dummy)
correct = dummy;
pattern = '[\./-]';

for i = 1:length(dummy)
    name = dummy{i};
    a = regexprep(name(1:end-1),pattern,'_');
    b = regexprep(name(end),pattern,'');
    correct{i} = [a,b];
end

