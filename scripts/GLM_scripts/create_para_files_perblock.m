%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create para condition timing files for
%%% the trifloc localizer data in order to run per-block regressions using
%%% freesurfer pipeline
%%% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/functions/');
ccc;

%% Initialize key variables/paths
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'}));
%subjCodes = {'MM', 'PP'};
N_subjs = length(subjCodes);

sections_per_block = 16;
base_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

% desired localizer condition order (sections_per_block=2):
% 1-4 = vA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 5-6 = vP (block 1, section 1; block 1, section 2)
% 7-10 = aA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 11-12 = aP (block 1, section 1; block 1, section 2)
% 13-16 = tA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 17-18 = tP (block 1, section 1; block 1, section 2)
% 19-20 = f (block 1, section 1; block 1, section 2)

% old:
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f

%% Loop through subjs and create para files for each run
for ss = 1:N_subjs
    subjCode = subjCodes{ss};
    subjdir = [base_dir subjCode '/localizer/'];
    files = dir(subjdir);
    N_runs = sum(contains({files.name}, '00'));
    disp(N_runs);
    assert(N_runs>=3, 'fewer than 3 runs');
    
    for rr = 1:N_runs
        old_parafile_path = [subjdir '00' num2str(rr) '/localizer_condition_timing.para'];
        old_para = readmatrix(old_parafile_path, 'FileType', 'text');
        new_para = old2new_para(old_para, sections_per_block);

        new_parafile_path = [subjdir '00' num2str(rr) '/localizer_condition_timing_perblock_16sections.para'];
        writematrix(new_para, new_parafile_path, 'filetype','text', 'delimiter','\t')
    end

end


function new_para = old2new_para(old_para, sections)
    
    N_conds = length(unique(old_para(:,2)));
    N_blocks = size(old_para,1);
    assert(N_blocks==12 | N_blocks==10, 'unexpected number of blocks (not 12 or 10)');

    % create condition lookup
    condition_convert_lookup = {};
    blocks_per_condition = arrayfun(@(x) sum(old_para(:,2)==x), 1:7, 'UniformOutput', false);
    for cc = 1:N_conds
        if cc == 1
            condition_convert_lookup{cc} = 1:(sections*blocks_per_condition{cc});
        else
            last_cond = condition_convert_lookup{cc-1}(end);
            condition_convert_lookup{cc} = last_cond+1:(last_cond+sections*blocks_per_condition{cc});
        end
    end

    % Make new para matrix
    new_para = nan(size(old_para,1)*sections, 3);
    for bb = 1:N_blocks
        blockdata = old_para(bb,:);
        for ss = 1:sections
            new_condition = condition_convert_lookup{blockdata(2)}(1);
            condition_convert_lookup{blockdata(2)} = condition_convert_lookup{blockdata(2)}(2:end); % remove condition from list after it has been used
            row_ind = (bb*sections)+(ss-sections);
            if ss == 1
                new_para(row_ind,:) = [blockdata(1) new_condition blockdata(3)/sections];
            else
                new_para(row_ind,:) = [new_para(row_ind-1,1)+new_para(row_ind-1,3) new_condition blockdata(3)/sections];
            end
        end
    end


end

