%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to count the number of runs available for
%%% each subj
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% 
subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N = length(subjCodes);
runs = nan(N,1);

for ss = 1:N
    subjCode = subjCodes{ss};
    runlist = readmatrix(['/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_runlistfile.txt']);
    runs(ss) = length(runlist)
end


