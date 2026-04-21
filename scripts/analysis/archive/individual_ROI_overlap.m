%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to assess the overlap between auditory,
%%% visual, and supramodal subject-level ROIs 
%%%
%%% Tom Possidente - January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROIs = {'aINS', 'aIPS', 'aMFG', 'ant_temporal', 'cIPS', 'CO', 'DO', 'inf_aINS', 'iPCS', 'LOT', 'midIFS', 'midINS', ...
    'ms_post_STSG', 'MT', 'parietal_opercular', 'pIPS', 'post_col_sulc', 'post_temporal', 'preSMA', 'sPCS',...
    'sup_aINS', 'tgPCS', 'VO', 'VOT'};
N_ROIs = length(ROIs);

hemis = {'lh', 'rh'};

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);

contrast_names = {'visual', 'auditory', 'supramodal'};

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific/';

%% Loop through subjs/hemis/ROIs and calculate overlap
dice_coefs = nan(N,2,N_ROIs,3); % subjs x hemis x ROIs x overlap combos

for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};
        for rr = 1:N_ROIs
            ROI_name = ROIs{rr};
            
            % Load ROI labels
            vis_ROI_fpath = [ROI_dir hemi '.' ROI_name '_visual_' subjCode '.label'];
            if isfile(vis_ROI_fpath)
                vis_ROI = readtable(vis_ROI_fpath, 'FileType','text');
            else
                vis_ROI = nan;
            end

            aud_ROI_fpath = [ROI_dir hemi '.' ROI_name '_auditory_' subjCode '.label'];
            if isfile(aud_ROI_fpath)
                aud_ROI = readtable(aud_ROI_fpath, 'FileType','text');
            else
                aud_ROI = nan;
            end
                        
            supra_ROI_fpath = [ROI_dir hemi '.' ROI_name '_supramodal_' subjCode '.label'];
            if isfile(supra_ROI_fpath)
                supra_ROI = readtable(supra_ROI_fpath, 'FileType','text');
            else
                supra_ROI = nan;
            end

            % Dice overlap between visual and auditory
            if istable(vis_ROI) && istable(aud_ROI)
                overlap_verts = sum(ismember(vis_ROI{:,1}, aud_ROI{:,1}));
                dice_coefs(ss, hh, rr, 1) = (2 * overlap_verts) / (height(vis_ROI)+height(aud_ROI)) ;
            end

            % Dice overlap between visual and supramodal
            if istable(vis_ROI) && istable(supra_ROI)
                overlap_verts = sum(ismember(vis_ROI{:,1}, supra_ROI{:,1}));
                dice_coefs(ss, hh, rr, 2) = (2 * overlap_verts) / (height(vis_ROI)+height(supra_ROI)) ;
            end

            % Dice overlap between auditory and supramodal
            if istable(aud_ROI) && istable(supra_ROI)
                overlap_verts = sum(ismember(aud_ROI{:,1}, supra_ROI{:,1}));
                dice_coefs(ss, hh, rr, 3) = (2 * overlap_verts) / (height(aud_ROI)+height(supra_ROI)) ;
            end
        end
    end
end


%% Plot dice coef results

for hh = 1:2
    data = squeeze(dice_coefs(:,hh,:,:));
    means = squeeze(mean(data, 1, 'omitnan'));
    Ns = squeeze(sum(~isnan(data),1));
    SEs = squeeze(std(data, 1, 'omitnan'))./sqrt(Ns);
    
    
    ngroups = size(means, 1);
    nbars = size(means, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    figure; hold on;
    bar(means);
    legend({'visual-auditory overlap', 'visual-supramodal overlap', 'auditory-supramodal overlap'});
    
    for bb = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*bb-1) * groupwidth / (2*nbars);
        errorbar(x, means(:,bb), SEs(:,bb), '.');
    end
    hold off
    xticks(1:N_ROIs);
    xticklabels(replace(ROIs, '_', ' '));
    xlabel('ROI');
    ylabel('Dice Coefficient');
    title(hemi);
end