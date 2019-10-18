function [machs,alts,sslips,mcenters] = getdata(blocks)
pattern1 = ' *\<MACH NO\> *= * (-*[\d]+.[\d]*)';
pattern2 = ' *\<ALTITUDE\> *= * (-*[\d]+.[\d]*)';
pattern3 = ' *\<SIDESLIP\> *= * (-*[\d]+.[\d]*)';
pattern4 = ' *\<MOMENT CENTER\> *= * (-*[\d]+.[\d]*)';

machs=[];
alts=[];
sslips=[];
mcenters=[];

for i = 1:length(blocks)
    block = blocks{i};
    mach = regexp(block,pattern1,'tokens');
    machs = [machs, str2double(mach{1})];
    sslip = regexp(block,pattern3,'tokens');
    sslips = [sslips, str2double(sslip{1})];
    alt = regexp(block,pattern2,'tokens');
    alts = [alts, str2double(alt{1})];
    mcenter = regexp(block,pattern4,'tokens');
    mcenters = [mcenters, str2double(mcenter{1})];
end