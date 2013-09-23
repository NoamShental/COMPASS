% Gerate figures for the paper:
% "Accurate Identification and Profiling of Microbial Communities using Massively Parallel Sequencing"
%

plot_identifiability=1;  % plot # of species identifiable as function of read length L. Figure 2 in paper
plot_error = 1;  % Plot reconstruction error as function of read length L. Figure 3 in paper

AssignGeneralConstants;

start_time = cputime;
switch machine
    case PC
        root_MCR_dir =  'C:\Users\orzuk\Google Drive\sync\Dropbox\bacterial_nextgen\';  % Change to your working directory
        %         input_dir = 'C:\Users\orzuk\Google Drive\sync\Dropbox\bacterial_nextgen\data\SimulationForRecomb';  % here are all seqeunces
        %         output_figs_dir = 'C:\\Users\orzuk\Google Drive\sync\Dropbox\bacterial_nextgen\figs\'; % save figures
        %         database_16s_file = '../../compressed_sensing/metagenomics/nextgen/data/bacteria_s16_data_uni.mat';
        %         database_16s_packed_file = '../../compressed_sensing/metagenomics/nextgen/data/s16_data_uni_packed.mat';
    case UNIX
        root_MCR_dir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/'; % Change to your working directory
        %         input_dir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/data/SimulationForRecomb';  % here are all seqeunces
        %         %        output_figs_dir = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/data/SimulationForRecomb/figs/'; % save figures
        %         database_16s_file = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/data/bacteria_s16_data_uni.mat';
        %         % database_16s_packed_file = '/seq/orzuk2/compressed_sensing/metagenomics/next_gen/data/s16_data_uni_packed.mat';
        %         database_16s_packed_file = '/seq/orzuk/identifiability/s16_data_uni_packed.mat'; % temp (due to space problems in /seq/orzuk2)
        %         output_figs_dir = '/seq/orzuk/identifiability/';  % temp (space problems)
end
input_dir = fullfile(root_MCR_dir, 'data/SimulationForRecomb'); % Take relative links
database_16s_file = fullfile(root_MCR_dir, 'bacteria_s16_data_uni.mat');
database_16s_packed_file = fullfile(root_MCR_dir, 's16_data_uni_packed.mat');
output_figs_dir = fullfile(root_MCR_dir, 'figs');

input_files = GetFileNames(fullfile(input_dir, 'For*.mat'), 1)
sparse_summary_file_name = fullfile(input_dir, 'all_runs_sparse_representation.mat');


num_files = length(input_files);
L = 100; % maximal read length


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Figure 3: Reconstruction Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(plot_error)
    S = load(database_16s_file, 'Header_uni', 'Sequence_packed', 'len_uni');
    num_reads_vec = zeros(num_files, 1); l2_error_vec = zeros(num_files, 1); mahalanobis_error_vec = zeros(num_files, 1);
    l2_error_std_vec = zeros(num_files, 1); mahalanobis_error_std_vec = zeros(num_files, 1);
    % l2_bound_vec = zeros(num_files, 1);
    SPRES = cell(num_files, 1); iters_vec = zeros(num_files, 1);
    iters = 100; %temp hard-coded
    l2_mat = zeros(num_files, iters);  mahalanobis_mat = zeros(num_files, iters); lambda_min_mat = zeros(num_files, iters);
    
    for i=1:num_files
        read_file = i
        if(~exist(sparse_summary_file_name, 'file'))
            RES = load(input_files{i});
            SPRES{i}.ovec = sparse(RES.ovec); SPRES{i}.rvec = sparse(RES.rvec); % make sparse
            iters = size(RES.ovec, 2);
        else
            if(i == 1)
                load(sparse_summary_file_name);
            end
            iters = size(SPRES{i}.ovec, 2); % allow different # of iterations for different runs
            iters_vec(i) = iters;
        end
        num_reads_vec(i) = abs(str2nums(input_files{i}));
        
        
        for j=1:iters
            l2_mat(i,j) = sqrt(sum((SPRES{i}.ovec(:,j) - SPRES{i}.rvec(:,j)).^2));
        end
        l2_error_vec(i) = mean(l2_mat(i,1:iters));
        l2_error_std_vec(i) = std(l2_mat(i,1:iters));
        
        
        for j=1:iters % Compute A^T A only on the two sets of species
            I_ovec = find(SPRES{i}.ovec(:,j)); I_rvec = find(SPRES{i}.rvec(:,j));
            
            
            build_matrix_j = j
            I = union(I_ovec, I_rvec); % get all indices
            [ReadBySpeciesMat kmers_packed] = ...
                BuildMixingMatrixFromSequences(L, S.Sequence_packed(I), S.len_uni(I), 32); % matlab word size (32 or 64??)
            
            A = ReadBySpeciesMat;
            for k=1:length(I)
                A(:,k) =  A(:,k) ./ (S.len_uni(I(k))-L+1);
            end
            D = A' * A;
            %            D = ReadBySpeciesMat' * ReadBySpeciesMat; % Matrix used for Mahalanobis
            mahalanobis_mat(i,j) = sqrt( (SPRES{i}.ovec(I,j) - SPRES{i}.rvec(I,j))' * D * (SPRES{i}.ovec(I,j) - SPRES{i}.rvec(I,j)) );
            lambda_min_mat(i,j) = min(eig(D)); % take minimal eigenvalue
        end % loop on iters
        mahalanobis_error_vec(i) = mean(mahalanobis_mat(i,1:iters));
        mahalanobis_error_std_vec(i) = std(mahalanobis_mat(i,1:iters));
    end
    
    delta = 1/2;
    dense_num_reads_vec = linspace(0, 10*max(num_reads_vec), 50000);
    
    mahalanobis_bound_vec = (2 + sqrt(1 / delta)) ./ sqrt(dense_num_reads_vec);
    lambda_min = 0.01;  % problem is we don't know lambda_min ...
    l2_bound_vec = mahalanobis_bound_vec ./ sqrt(lambda_min);
    
    if(~exist(sparse_summary_file_name, 'file'))
        save(sparse_summary_file_name, 'SPRES', 'iters', 'num_reads_vec', 'iters_vec', ...
            'l2_mat', 'mahalanobis_mat', 'lambda_min_mat', ...
            'l2_error_vec', 'l2_error_std_vec', 'l2_bound_vec', ...
            'mahalanobis_error_vec', 'mahalanobis_error_std_vec', 'mahalanobis_bound_vec');
    end
    
    mahalanobis_bound_vec = mahalanobis_bound_vec .* sqrt(mean(S.len_uni));
    mahalanobis_error_vec = mahalanobis_error_vec .* sqrt(mean(S.len_uni));
    
    
    for log_x_flag = 0:1 % plot in linear and log scales
        for log_y_flag = 0:1
            fig = figure;
            switch log_x_flag
                case 0
                    log_x_str = '';
                    switch log_y_flag
                        case 0
                            log_y_str = '';
                            errorbar(num_reads_vec, l2_error_vec, l2_error_std_vec, '*', 'linewidth', 2); hold on;
                            plot(dense_num_reads_vec, l2_bound_vec, 'linewidth', 2); hold on;
                            errorbar(num_reads_vec, mahalanobis_error_vec, mahalanobis_error_std_vec, '*r', 'linewidth', 2); hold on;
                            plot(dense_num_reads_vec, mahalanobis_bound_vec, 'r', 'linewidth', 2); hold on;
                        case 1 % semilogy
                            log_y_str = '_log_y';
                            errorbar(num_reads_vec, l2_error_vec, l2_error_std_vec, '*', 'linewidth', 2); errorbarlogy; hold on;
                            semilogy(dense_num_reads_vec, l2_bound_vec, 'linewidth', 2); hold on;
                            errorbar(num_reads_vec, mahalanobis_error_vec, ...
                                min(mahalanobis_error_vec-0.001, mahalanobis_error_std_vec), ...
                                mahalanobis_error_std_vec,'*r', 'linewidth', 2); hold on;
                            semilogy(dense_num_reads_vec, mahalanobis_bound_vec, 'r', 'linewidth', 2); hold on;
                            
                    end % switch on y flag
                case 1
                    log_x_str = '_log_x';
                    switch log_y_flag
                        case 0 % semilogx
                            log_y_str = '';
                            errorbar(num_reads_vec, l2_error_vec, l2_error_std_vec, '*', 'linewidth', 2); errorbarlogx; hold on;
                            semilogx(dense_num_reads_vec, l2_bound_vec, 'linewidth', 2); hold on;
                            errorbar(num_reads_vec, mahalanobis_error_vec, mahalanobis_error_std_vec, '*r', 'linewidth', 2); errorbarlogx; hold on;
                            semilogx(dense_num_reads_vec, mahalanobis_bound_vec, 'r', 'linewidth', 2); hold on;
                        case 1 % loglog
                            log_y_str = '_log_y';
                            errorbar(num_reads_vec, l2_error_vec, l2_error_std_vec, '*', 'linewidth', 2);
                            ax = get(fig,'CurrentAxes'); hold on;
                            set(ax,'XScale','log','YScale','log')
                            loglog(dense_num_reads_vec, l2_bound_vec, 'linewidth', 2); hold on;
                            errorbar(num_reads_vec, mahalanobis_error_vec, ...
                                min(mahalanobis_error_vec-0.00000001, mahalanobis_error_std_vec), ...
                                mahalanobis_error_std_vec, '*r', 'linewidth', 2); hold on;
                            loglog(dense_num_reads_vec, mahalanobis_bound_vec, 'r', 'linewidth', 2); hold on;
                            y_lim = ylim; ylim([y_lim(1) 1]); xlim([5000 5*10^6]);
                    end % switch log_y_flag
            end % switch log_x_flag
            xlabel('# Reads', 'fontsize', 16, 'fontweight', 'bold'); ylabel('Error', 'fontsize', 16, 'fontweight', 'bold');
            legend({'L_2 (simulations)', 'L_2 (upper-bound)', 'MA (simulations)', 'MA (upper-bound)'}, 'fontsize', 16, 'fontweight', 'bold');
            set(gca, 'fontsize', 16);   legend boxoff;
            
            %%%        my_saveas(gcf, fullfile(input_dir, 'figs', ['L_2_error_as_function_of_num_reads' log_x_str log_y_str]), {'pdf', 'epsc'});
            
        end % loop on y flag
    end % loop on x flag
end % if plot error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Figure 2: Identifiability as function of read length %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(plot_identifiability)
    %    for iter = [1 2 3 4 6] % loop on finished iterations
    if(~exist('iter', 'var'))
        iter=1;
    end
    if(~exist('k_species', 'var'))
        k_species = 900;
    end
    if(~exist('max_L', 'var'))
        max_L=30; % L_vec = 1:100;
    end
    if(~exist('min_L', 'var'))
        min_L= 10; % don't go to too low read length (causes problems)
    end
    if(~exist('rand_flag', 'var'))
        rand_flag = 1;  % take random species
    end
    if(rand_flag)
        rand_str = '_rand';
    else
        rand_str = '';
    end
    output_data_file = fullfile(output_figs_dir, ...
        ['identifiability_vs_read_length_' num2str(k_species) '_species_iter_' num2str(iter) ...
        '_L_' num2str(min_L) '_' num2str(max_L) rand_str '.mat'])
    if(machine == UNIX) % run and compute
        S = load(database_16s_packed_file, 'Header_uni', 'Sequence_packed', 'len_uni', 'bad_seq_inds');
        num_species = length(S.len_uni);
        rng('shuffle'); % make sure that the rng is changed in each run
        S.good_seq_inds = setdiff(1:num_species, S.bad_seq_inds)
        if(rand_flag)
            I = randperm(num_species); I = I(1:k_species); % pick at random 10000 species
        else
            I=S.good_seq_inds(1:k_species);  % I=1:k_species; % Set fixed!! (first k) %
        end
        num_identifiable = zeros(max_L, 1); rand_num_identifiable = zeros(max_L,1);
        identifiable_inds = cell(max_L, 1); rand_identifiable_inds = cell(max_L, 1);
        rand_seqs = cell(k_species,1); rand_seqs_packed = cell(k_species,1);
        for i=1:k_species
            rand_seqs{i} = randi(4, 1, S.len_uni(I(i)));
            rand_seqs_packed{i} = pack_seqs(rand_seqs{i}, matlab_word_size);
        end
        for L=min_L:max_L % :-1:min_L
            check_identifiability_L = L
            [identifiable_inds{L} A kmers_packed] = get_identifiable_species(S.Sequence_packed(I), S.len_uni(I), L); % for each read length - find which species are identifiable
            num_identifiable(L) = sum(identifiable_inds{L}); % we only care about the number of identifiable species
            
            if(rand_flag) % NEW! get also random sequences
                [rand_identifiable_inds{L} rand_A rand_kmers_packed] = get_identifiable_species(rand_seqs_packed, S.len_uni(I), L); % for each read length - find which species are identifiable
                rand_num_identifiable(L) = sum(rand_identifiable_inds{L}); % we only care about the number of identifiable species
            end
        end
        save(output_data_file, 'min_L', 'max_L', 'k_species', 'rand_flag', ...
            'num_identifiable', 'I', 'identifiable_inds', 'rand_num_identifiable', 'rand_identifiable_inds', ...
            'A', 'kmers_packed'); % save figure data
    end
    if(machine == PC) % here machine is PC
        output_figs_dir = 'identifiability\';
        num_identifiable_vec=zeros(100,1); rand_num_identifiable_vec=zeros(100,1);
        for min_L = [5:9 10 15:10:95] % should be 5 to 95 % first load data
            if(min_L < 10)
                max_L=min_L;
            else
                if(min_L==10)
                    max_L=min_L+4;
                else
                    max_L=min_L+9;
                end
            end
            output_data_file = fullfile(output_figs_dir, ...
                ['identifiability_vs_read_length_' num2str(k_species) '_species_iter_' num2str(iter) ...
                '_L_' num2str(min_L) '_' num2str(max_L) rand_str '.mat'])
            UUU{min_L} = load(output_data_file, 'num_identifiable'); % load figure data
            try
                UUU_rand{min_L} = load([remove_suffix_from_file_name(output_data_file) '_rand.mat'], 'rand_num_identifiable');
            catch
                UUU_rand{min_L}.rand_num_identifiable(min_L:max_L) = max(rand_num_identifiable_vec); % -10; % fill with -1's to make sure we know we didn't load
            end
            num_identifiable_vec(min_L:max_L) = UUU{min_L}.num_identifiable(min_L:max_L);
            rand_num_identifiable_vec(min_L:max_L) = UUU_rand{min_L}.rand_num_identifiable(min_L:max_L);
        end
        
        close all;
        big_fig = figure(1); plot(1:max_L, num_identifiable_vec ./ k_species, 'r', 'linewidth', 2); hold on;
        plot(1:max_L, rand_num_identifiable_vec ./ k_species, 'b', 'linewidth', 2);
        xlabel('Read Length'); ylabel('Frac. Identifiable Species'); xlim([1 100]); ylim([0 1.005]);
        %        title(['Identifiability for ' num2str(k_species) ' species']);
        inset_fig = figure(2); plot(1:max_L, num_identifiable_vec ./ k_species, 'r', 'linewidth', 2); hold on;
        plot(1:max_L, rand_num_identifiable_vec ./ k_species, 'b', 'linewidth', 2);
        %        xlabel('Read Length'); ylabel('Frac. Identifiable Species');
        xlim([7 100]); ylim([0.96 1.001]);
        
        [h_m h_i]=inset(big_fig,inset_fig,0.7);
        %        set(h_i,'xtick',2.35:.025:2.45,'xlim',[2.35,2.45])
        
        legend(h_m, {'16s', 'rand.'}, 4); legend(h_m, 'boxoff');        
        my_saveas(gcf, fullfile(output_figs_dir, 'identifiability_vs_read_length'), {'epsc', 'pdf'}); % save figure 
    end % if machine
    %    end % loop on iterations
end % if plot identifiability




% Run jobs:
for rand_flag = 1:1 % 0:0 % 1
    for min_L=5:20 % 5:10:95
        max_L=min_L; % +15;
        for iter=1:(1+3*rand_flag) % need many iterations only for random
            for k_species_vec = 10000 % [5000 10000 20000 50000] % 500 1000
                k_species = k_species_vec
                if(rand_flag)
                    rand_str = '_rand';
                else
                    rand_str = '';
                end
                
                if(k_species_vec>10000)
                    queue_str = 'priority';
                else
                    queue_str = 'hour';
                end
                job_str = ['rand_flag=' num2str(rand_flag) '; iter=' num2str(iter) '; k_species=' num2str(k_species_vec) '; min_L=' num2str(min_L) '; max_L=' num2str(max_L) '; make_fig_error_as_function_of_read_length;']
                SubmitMatlabJobToFarm(job_str, ...
                    ['out/identifiability_' num2str(k_species) '_species_iter' num2str(iter) ...
                    '_L_' num2str(min_L) '_' num2str(max_L) rand_str '.out'], queue_str, [], [], 12);
            end
        end
    end
end

make_figs_time = cputime - start_time
