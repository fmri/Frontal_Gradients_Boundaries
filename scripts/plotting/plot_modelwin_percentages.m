%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to visualize individual subj model win
%%% percentages (gradient, boundary, inconclusive).
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROIs = {'preSMA', 'PMv', 'aINS', 'PMd', 'aIPS', 'midDLPFC'};
hemis = {'lh', 'rh'};
modalities = {'visual', 'auditory'};

data = table({'visual'; 'auditory'; ''; 'visual'; 'auditory'; ''; 'visual'; 'auditory'; ''; 'visual'; 'auditory'; ''; 'visual'; 'auditory'; ''; 'visual'; 'auditory'},....
             {'preSMA'; 'preSMA'; ''; 'PMv'; 'PMv'; ''; 'aINS'; 'aINS'; ''; 'PMd'; 'PMd'; ''; 'aIPS'; 'aIPS'; ''; 'midDLPFC'; 'midDLPFC'},...
             {'combined'; 'combined'; ''; 'combined'; 'combined'; ''; 'left'; 'combined'; ''; 'left'; 'combined'; ''; 'none'; 'right'; ''; 'left'; 'left'},...
             [77; 69; 0; 83; 67; 0; 63; 63; 0; 61; 27; 0; 0; 53; 0; 33; 25],...
             [11; 17; 0; 6; 7; 0; 25; 14; 0; 6; 15; 0; 0; 20; 0; 53; 50],...
             [12; 14; 0; 11; 26; 0; 12; 23; 0; 33; 58; 0; 0; 27; 0; 14; 25], ...
             'VariableNames',{'modality', 'ROI', 'hemisphere', 'gradient', 'boundary', 'inconclusive'});

%% Plot
figure; 
bar(data{:,4:6});
xticks([1,2,4,5,7,8,10,11,13,14,16,17]);
xticklabels({'preSMA vis C', 'preSMA aud C', 'PMv vis C', 'PMv aud C',...
            'aINS vis L', 'aINS aud C', 'PMd vis L', 'PMd aud C',...
            'aIPS vis', 'aIPS aud R', 'midDLPFC vis L', 'midDLPFC aud L'});
legend({'Gradient', 'Boundary', 'Inconclusive'});
ylabel('Percent Win');
ylim([0,100]);
fontsize(16,"points")