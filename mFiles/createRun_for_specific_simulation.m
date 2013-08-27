% the same as example_of_a_single_simulation - but receives parameters through the struct auxDataIn

function [uniqueReads,uniqueReads_length,correctWeight]=createRun_for_specific_simulation(auxDataIn)

% read errrors 
ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.005; 
ErrorStruct.final_error = 0.03; 

p_one_nuc_error_mat = [ -1   0.3  0.22  0.18
                    0.5    -1  0.22   0.6
                    0.35  0.15    -1  0.22
                    0.15  0.55  0.56   -1]; 

ErrorStruct.p_one_nuc_error_mat = p_one_nuc_error_mat;


% create the mixture
ind_bac_in_mix=randi([1,auxDataIn.numBacteriaInDatabase],auxDataIn.numberOfBacteriaInMix,1);
w=1./([1:auxDataIn.numberOfBacteriaInMix].^auxDataIn.powerLaw);
w = w./sum(w);

correctWeight = zeros(1,auxDataIn.numBacteriaInDatabase);
correctWeight(ind_bac_in_mix) = w;


% create reads
auxData = struct;
auxData.correctWeight = correctWeight; % vector 
auxData.Nreads = auxDataIn.Nreads;
auxData.readLength = auxDataIn.readLength;
auxData.seed_in = 1;% seed for reads creation
auxData.basicSeqNameDir = [auxDataIn.userDir,'database/full16S/datNoNonACGT/packed64/'];
auxData.basicSeqKey= [auxDataIn.userDir,'database/full16S/datNoNonACGT/keyNoNonACGT'];
auxData.addNoiseFlag = auxDataIn.addNoiseFlag; % 1 for adding noise, 0= no noise
auxData.ErrorStruct = ErrorStruct;

% create reads in a packed format
[uniqueReads,uniqueReads_length] = createReads_package(auxData);

