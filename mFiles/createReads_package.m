% simulate reads according a set of parameters in auxData.
% the output are a list of unique reads and the number each of the appears. % Each read is packed as a 64bit word
% For an example see example_of_a_single_simulation.m


function [uniqueReads,uniqueReads_length]=createReads_package(auxData)

% randomize seed
if ~isfield(auxData,'seed_in')
  outseed = sum(100*clock);
else
  outseed = auxData.seed_in;  
end
rand('seed',outseed);
%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the parametes
ind_bac_in_mix = find(auxData.correctWeight); % indices of bacteria in the mixture
freq = auxData.correctWeight(ind_bac_in_mix); % frequency of bacteria
numReadsBacteria = round(freq/sum(freq)*auxData.Nreads); % number of reads
Nbac_in_mixture = length(ind_bac_in_mix); % number of bacteria in the mixture

% load sequences of bacteria in the mixture.
load(auxData.basicSeqKey,'positionInPart','len_uni');
len_mix = len_uni(ind_bac_in_mix); clear len_uni
Sequence_mix = cell(length(ind_bac_in_mix),1);
for i=1:length(ind_bac_in_mix)
  clear seq_*
  load([auxData.basicSeqNameDir,'seq_part_',num2str(positionInPart(ind_bac_in_mix(i)))],['seq_',num2str(ind_bac_in_mix(i))]);
  ww = ['Sequence_mix{i} = seq_',num2str(ind_bac_in_mix(i)),'{1};'];
  eval(ww);
end
clear seq_* positionInPart 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



% create the reads 
currSeq = extract_sub_kmers(Sequence_mix{1}, len_mix(1), auxData.readLength, 0,0,1,1); % simulate one read, a stupid way of calculating the size of the vector in its packed form.

simulated_reads = zeros(auxData.Nreads,length(currSeq),'uint64');clear currSeq
k = 1;
%keyboard
for i=1:Nbac_in_mixture
  bac_read_startvec=randi([1,len_mix(i)-auxData.readLength+1],numReadsBacteria(i),1); % select the reads' random starting points for the i'th bacterium 
  for j=1:numReadsBacteria(i) 
    simulated_reads(k+j-1,:) = extract_sub_kmers(Sequence_mix{i}, len_mix(i), auxData.readLength, 0,0,1,bac_read_startvec(j)); % create the read
  end
  k=k+numReadsBacteria(i);  
end
simulated_reads(k:end,:) = [];

% add noise if needed
if auxData.addNoiseFlag
  read_error_substitution_table = GenerateSubstitutionErrorTable(auxData.readLength, auxData.ErrorStruct);
  
  if size(simulated_reads,1)>10^5 % if more than 10^5 reads - add noise for each group of 10^5 separately - for memory purposes
    partRed = 1:10^5:size(simulated_reads,1);
    partRed(end+1) = size(simulated_reads,1)+1;
    
    for ind_part_red=1:length(partRed)-1
      curPartRed = partRed(ind_part_red):partRed(ind_part_red+1)-1;
      simulated_reads(curPartRed,:) = add_noise_to_kmers(auxData.readLength, size(curPartRed,2), ...
        simulated_reads(curPartRed,:), reshape(read_error_substitution_table, auxData.readLength, 4*4));
    end
  else % for less than 100000
    partRed = 1;
    partRed(end+1) = size(simulated_reads,1)+1;
    for ind_part_red=1:length(partRed)-1
      curPartRed = partRed(ind_part_red):partRed(ind_part_red+1)-1;
      simulated_reads(curPartRed,:) = add_noise_to_kmers(auxData.readLength, size(curPartRed,2), ...
        simulated_reads(curPartRed,:), reshape(read_error_substitution_table, auxData.readLength, 4*4));
    end
    
  end
  
  
  %keyboard
  disp(['added noised according to: model: ',auxData.ErrorStruct.error_model,' baseline_error: ',num2str(auxData.ErrorStruct.baseline_error),' final_error: ',num2str(auxData.ErrorStruct.final_error)]);

end

% end creating reads
%%%%%%%%%%%%%%%%%%%%%%%

%keyboard

% create a list of unique reads and pack them in 64bit words in case they are not already packed as on this case
[uniqueReads,uniqueReads_length] = createUniqueReadsAndPack(simulated_reads,auxData.readLength);



