% NormalizeReads
% input : 
% filename - the k-mers (after using a running window of 90) and their
% average mapping position (from NOAM????? program)
% globalFilterWinSize - the window size for the mean filter used to
% calculate the global noise (we used 10)
% globalFilterStartPos - the first position to be mean filtered (so we
% don't insert the initial step to the filter) (we used 40)
%
% output:
% uniqueReads - the reads themselves - the same as the input file.
% uniqueReads_length - the normalized counts per k-mer
% uniqueReads_length_before_normalization - the original counts from the input file - used in the postprossing step
% globalProfileForward - the global,mean-filtered noise profile for forward reads
% globalProfileReverse - the global,mean-filtered noise profile for reverse reads

function [uniqueReads,uniqueReads_length,uniqueReads_length_before_normalization,globalProfileForward,globalProfileReverse]=normalizeExperimentalReads(filename,globalFilterWinSize,globalFilterStartPos)

%% Calculate noise profiles
avgDat=load(filename);

% prepare the number of reads per position buffers (forward and reverse)
MAX16SSIZE=1600;
valFor=zeros(MAX16SSIZE,1);
valRev=zeros(MAX16SSIZE,1);
valRevRev=zeros(MAX16SSIZE,1);
forRevRatio=zeros(MAX16SSIZE,1);
pos=1:MAX16SSIZE;

% go over all reads 
% for each read use the mapping to the avg. position to create a histogram
for a=1:length(avgDat.uniqueReads)
    ccell=avgDat.shift_and_POS_data{a};
    for b=1:length(ccell.shift)
        cpos=round(ccell.POS(b)+ccell.shift(b)-1);
        cposr=round(ccell.REV_POS(b)+11-ccell.shift(b));

        valFor(cpos)=valFor(cpos)+ccell.length_forward(b);
        valRev(cposr)=valRev(cposr)+ccell.length_reverse(b);
        % and valRevRev for the ratio
        valRevRev(cpos)=valRevRev(cpos)+ccell.length_reverse(b);
    end
end

% calculate the forward/reverse ratio
for a=1:length(valFor)
    forRevRatio(a)=valFor(a)/(valFor(a)+valRevRev(a));
end

%% Smooth using mean filter length WINSIZE
WINSIZE=globalFilterWinSize;
% STARTPOS - to not include the initial spike in the normalization
STARTPOS=globalFilterStartPos;

globalProfileForward=valFor;
for a=STARTPOS:length(valFor)
    globalProfileForward(a)=mean(valFor(max((a-WINSIZE),STARTPOS):min((a+WINSIZE),length(globalProfileForward))));
end
globalProfileReverse=valRev;
for a=STARTPOS:length(valFor)
    globalProfileReverse(a)=mean(valRev(max((a-WINSIZE),STARTPOS):min((a+WINSIZE),length(globalProfileReverse))));
end

%% Normalize reads
uniqueReads_length=zeros(size(avgDat.uniqueReads_length_forward));
for a=1:size(uniqueReads_length,1)
    ccell=avgDat.shift_and_POS_data{a};
    cavg=0;
    cavgnum=0;
    for b=1:length(ccell.shift)
        cpos=round(ccell.POS(b)+ccell.shift(b)-1);
        cposr=round(ccell.REV_POS(b)+11-ccell.shift(b));

        
        if (globalProfileForward(cpos)>0)
            uniqueReads_length(a)=uniqueReads_length(a)+(ccell.length_forward(b)/globalProfileForward(cpos))*forRevRatio(cpos);
        end

        cavg=cavg+cpos;
        cavgnum=cavgnum+1;

        
        if (globalProfileReverse(cposr)>0)
            uniqueReads_length(a)=uniqueReads_length(a)+(ccell.length_reverse(b)/globalProfileReverse(cposr))*(1-forRevRatio(cpos));
        end
        
    end
end

uniqueReads_length = ceil(uniqueReads_length*10^5); % multiply by 10^5 to have integers. 
uniqueReads = avgDat.uniqueReads;
uniqueReads_length_before_normalization = avgDat.uniqueReads_length_forward+avgDat.uniqueReads_length_reverse; % the original nmuber of times each
                                                                                                               % read appears

