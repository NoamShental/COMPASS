% an example of loading reads from a fasta file and converting them to a COMPASS format

userDir = '~/CS/mFilesBAC/COMPASS/'; % ENTER THE PATH TO THE COMPASS DIRECTORY

% load the fasta file
% the file should include reads of constant length - Hence first filter the reads and trim them to constant length
fastaFileName = [userDir,'experimentalReads/exampleFastaFile.fa'];
[header,seq] = fastaread(fastaFileName);
seq_matrixForm = cell2mat(seq');

readLength = 90; % the read length

% COMPASS reads 
[uniqueReads,uniqueReads_length]=createUniqueReadsAndPack(seq_matrixForm,readLength);

% uniqueReads is a list of unique reads. Each read is saved as a 64bit word.
% uniqueReads_length is the number of times each such read appears

