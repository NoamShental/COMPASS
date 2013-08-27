% calculate mothur pairwise sequence distances between original fastq and
% recreated fasta (for emirge comparison)
% input:
% dirname - the direcory where the data is
% fastqname - the sequences from the original mixture in fastq format
% fastaname - the recreated emirge fasta file
% outfilename - the file to write distances to

function [distmat,numfq,numfa]=CalcFastaDist(dirname,fastqname,fastaname,outfilename)

pwd1 = pwd;
%cd(dirname);
%keyboard
% join the fastq and fasta
[numfa,numfq]=JoinFastqFasta(dirname,fastqname,fastaname,[outfilename '.fa']);
disp('Original');
disp(numfq);
disp('recreated');
disp(numfa);


% run the mothur for multiple sequence alignment and distance matrix
% calculation
disp('starting mothur');
%keyboard

cd(dirname);

basedir='~';

w = ['!',basedir,'/CS/BAC/mothur/mothur "#align.seqs(candidate=',...
     outfilename,'.fa'...
     ', template=',basedir,'/CS/BAC/core_set_aligned.fasta.imputed); dist.seqs(fasta=',...
     outfilename,'.align, output=square)"'];
%w = ['!',auxData.userDir,'/CS/BAC/mothur/mothur "#align.seqs(candidate=',...
%     filename,'.fa'...
%     ', template=',auxData.userDir,'/CS/BAC/core_set_aligned.fasta.imputed); dist.seqs(fasta=',...
%     filename,'.align, output=square)"'];
eval(w)

disp('Finished calculating distances');

% load the mothur output matrix and skip the first line
%distmat=dlmread([auxData.currDir,'/',filename,'.square.dist'],'\t',0,1);
distmat=dlmread([outfilename,'.square.dist'],'\t',0,1);
distmat=distmat(2:end,1:end-1);

% save the output
%save([fileName,'_mothurRes'],'distmat','allSeqs','origFreqVec','recFreqVec');
save([outfilename,'_mothurRes'],'distmat','fastqname','fastaname','numfa','numfq','dirname');


%unix(['cd ',auxData.currDir,';rm -f ',filename,'.align.report ',tmpfilename,'.align ',tmpfilename,'.square.dist mothur.*.logfile'])
unix(['rm -f ',outfilename,'.align.report ',outfilename,'.align ',outfilename,'.square.dist mothur.*.logfile']);
cd(pwd1);
disp('Minimal distances:');
disp(min(distmat(1:numfq,numfq+1:end)'));
disp('num identical:');
disp(sum(min(distmat(1:numfq,numfq+1:end)')==0));
disp('num close (2%):');
disp(sum(min(distmat(1:numfq,numfq+1:end)')<=0.02));
