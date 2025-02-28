%% Behavior Analyse from MedPC *.csv files

% The code relies on the following adapted function:
% --> ZeroOnesCount.m
%     https://www.mathworks.com/matlabcentral/answers/16356-finding-sequences-of-1-s-values 


% Variable -> data.bevaviour or data.bevaviour_bin (Averaging time epochs at desired sample rate)
% - first  cell row
%   .full record.

% - second cell row - from full session
%   . Row 1: index values for each event 
%   . Row 2: time in samples for each event
%   . Row 3: percentage for each event

% - third cell row
%   . Column 1: Baseline
%   . Column 2: CS-Trial | Inter trial periods (ITI-epoch)

% - fourth cell line
%   . Imobillity lower threshold. Binary: 1 if <= threshold & 0 if > threshold

% - fifth cell line
%   . Row 1: index values for each event
%   . Row 2: time in samples for each event
%   . Row 3: percentage for each event

% - sixth cell line
%   . Total time in sec

% - seventh cell line
%   . Percentage %

    
% by Flavio Mourao.
% email: mouraofg@illinois.edu
% Maren Lab - Beckman Institute for Advanced Science and Technology
% TUniversity of Illinois Urbana-Champaign

% Started in:  12/2023
% Last update: 12/2023

%% Load and convert table to array
data.behavior_MedPc =[];
%data.behavior_MedPc{1,1} = table2array(SSC6freezingextinction)';
data.behavior_MedPc{1,1} = table2array(Extinctionday1movementdata)';
data.behavior_MedPc{1,1} = Behavior_raw_extinction';

%% Interpolate data from MedPc to Miniscope sample size

parameters.sampling_rate_MedPC = 5;
% 
% factor = parameters.sampling_rate/parameters.sampling_rate_MedPC;
% data.behavior =[];
% data.behavior{1,1} = interp(data.behavior_MedPc{1,1},factor);
% 
% clear('SSC6freezingextinction','factor')

%% Data from behavior considering MedPC sample rate.
% Trial epochs - baseline and CS sound period
% Original sample rate

% Freezing settings
parameters.thr_1 = 5; % lower threshold in percentage.
parameters.thr_2 = 1; % higher threshold in sec to consider freezing.

% Full session Freezing

for ii = 1:size(data.behavior_MedPc{1,1},1)

    % idx for freezing epochs
    [data.behavior_MedPc{ii+1,1}(1,:), data.behavior_MedPc{ii+1,1}(2,:), ~] = ZeroOnesCount(data.behavior_MedPc{1,1}(ii,:) <= parameters.thr_1); % first row = start events / second = time in samples
    idx_to_remove = data.behavior_MedPc{ii+1, 1}(2,:) < parameters.thr_2 * parameters.sampling_rate_MedPC; % Exclude values based on timing threshold
    data.behavior_MedPc{ii+1, 1}(:,idx_to_remove) = [];
    data.behavior_MedPc{ii+1, 1}(3,:) = data.behavior_MedPc{ii+1, 1}(2,:)./parameters.sampling_rate_MedPC; % thirt row = time in second


    % idx for non-freezing epochs
    % First lets find all freezing idx
    f_idx(:,1) = data.behavior_MedPc{ii+1, 1}(1,:);
    f_idx(:,2) = data.behavior_MedPc{ii+1, 1}(1,:) + data.behavior_MedPc{ii+1, 1}(2,:);

    %Then, lets find the intervals who represents non-freezing
    data.behavior_MedPc{ii+1,2}(:,1) = [1 f_idx(1,1)-1]; % from beggining until the first freezing event
    data.behavior_MedPc{ii+1,2}(1,2:length(f_idx)) = f_idx(1:end-1,2)+1;                      % first row = start events
    data.behavior_MedPc{ii+1,2}(2,2:length(f_idx)) = (f_idx(2:end,1)-1) - f_idx(1:end-1,2)-1; % second = time in samples
    data.behavior_MedPc{ii+1,2}(3,:) = data.behavior_MedPc{ii+1, 2}(2,:)./parameters.sampling_rate_MedPC;   % thirt row = time in second

    idx_to_remove_1 = data.behavior_MedPc{ii+1, 2}(2,:) < parameters.thr_2 * parameters.sampling_rate_MedPC; % Exclude values based on timing threshold
    data.behavior_MedPc{ii+1, 2}(:,idx_to_remove_1) = [];

    clear('f_idx','idx_to_remove','idx_to_remove_1')

end


%%

initiation_samples  = zeros(max(max(cell2mat(cellfun(@length,data.behavior_MedPc(2:end,:),'UniformOutput',false)))),23);
termination_samples = zeros(max(max(cell2mat(cellfun(@length,data.behavior_MedPc(2:end,:),'UniformOutput',false)))),23);

for ii = 2:size(data.behavior_MedPc,1)
    
    idx_1 = data.behavior_MedPc{ii,1}(1,:)';
    initiation_samples(1:length(idx_1),ii-1) = idx_1;

    idx_2 = (data.behavior_MedPc{ii,1}(1,:) + data.behavior_MedPc{ii,1}(2,:))';
    termination_samples(1:length(idx_2),ii-1) = idx_2;

end

initiation_sec  = zeros(max(max(cell2mat(cellfun(@length,data.behavior_MedPc(2:end,:),'UniformOutput',false)))),23);
termination_sec = zeros(max(max(cell2mat(cellfun(@length,data.behavior_MedPc(2:end,:),'UniformOutput',false)))),23);


for ii = 2:size(data.behavior_MedPc,1)
    
    idx_1 = (data.behavior_MedPc{ii,1}(1,:)./ parameters.sampling_rate_MedPC)';
    initiation_sec(1:length(idx_1),ii-1) = idx_1;

    idx_2 = ((data.behavior_MedPc{ii,1}(1,:) + data.behavior_MedPc{ii,1}(2,:)) ./ parameters.sampling_rate_MedPC)';
    termination_sec(1:length(idx_2),ii-1) = idx_2;

end

%% test plot

figure
set(gcf,'color','white')

t = {'MT1','MT2','MT3','MT4','MT5','MT6','MT7','MT8','MT9','MT10','BUZ1','BUZ3','BUZ4','Rerec','RErec1','RErec3','RErec4',...
    'RErec5','RErec6','RErec7','RErec8','RErec10','RErec11'};

for ii = 1:size(data.behavior_MedPc{1,1},1)

    subplot(5,5,ii)
    
    plot(data.behavior_MedPc{1, 1}(ii,:),'color',[.6 .6 .6])
    hold on
    fill([initiation_samples(:,ii) termination_samples(:,ii) termination_samples(:,ii) initiation_samples(:,ii)], [0 0 100 100], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    title(t{1,ii})

    yline(5,'r','linewidth',2)
    xlim([1 10000])
    box off

end


%% Trial epochs - Extinction, Retrieval - baseline, CS and ITI epochs

% baseline
data.behavior_MedPc{3,1}  = data.behavior_MedPc{1,1}(1,1: round(data.events{3, 1}(1,1) .* parameters.sampling_rate_MedPC) -1);

epochs = [];
% CS sound and ITI period
for ii = 1:size(data.events{2, 1},1)

        epochs{1,ii}  = data.behavior_MedPc{1,1}(1, round(data.events{3, 1}(ii,1).*parameters.sampling_rate_MedPC) : round(data.events{3, 1}(ii,2).*parameters.sampling_rate_MedPC)); % CS-Trials
        epochs{2,ii}  = data.behavior_MedPc{1,1}(1, round(data.events{3, 2}(ii,1).*parameters.sampling_rate_MedPC) : round(data.events{3, 2}(ii,2).*parameters.sampling_rate_MedPC)); % ITI-Trials

end

% Reshaped CS and ITI trials in correct order and add to data behavior epochs
CS_ITI_Trials = reshape(epochs,1,[]);
data.behavior_MedPc(3,2:length(CS_ITI_Trials)+1) = CS_ITI_Trials;


% CS-Trials and ITI Freezing
for ii = 1:size(data.behavior_MedPc,2)
    
    data.behavior_MedPc{4,ii} = data.behavior_MedPc{3,ii} <= parameters.thr_1;

    [data.behavior_MedPc{5,ii}(1,:), data.behavior_MedPc{5,ii}(2,:), ~] = ZeroOnesCount(data.behavior_MedPc{4,ii});

    idx_to_remove = data.behavior_MedPc{5, ii}(2,:) < parameters.thr_2 * parameters.sampling_rate_MedPC;
    data.behavior_MedPc{5, ii}(:,idx_to_remove) = [];
    data.behavior_MedPc{5, ii}(3,:) = data.behavior_MedPc{5, ii}(2,:)./parameters.sampling_rate_MedPC;

    idx_to_remove = [];

    data.behavior_MedPc{6, ii} = sum(data.behavior_MedPc{5, ii}(3,:));
    data.behavior_MedPc{7, ii} = (data.behavior_MedPc{6, ii}.*100) ./ (length(data.behavior_MedPc{3, ii}) ./ (parameters.sampling_rate_MedPC));
end

clear('ii','epochs','CS_ITI_Trials','parameters.thr_1','parameters.thr_2','idx_to_remove')

%% Graph based on Maren`s lab
% Data freezing related to extinction S. Following the laboratory method, averages were calculated for every 5 time windows.

% Extinction
% Choose data:
data_2_use = data.behavior_MedPc(:,1:end);

data.behavior_MedPc_CS_epochs = zeros(1,(round(length(data_2_use)/5))/2 +1); 
data.behavior_MedPc_CS_epochs(1,1) = cell2mat(data_2_use(7, 1)); % baseline
data.behavior_MedPc_CS_epochs(1,2:end) = mean(reshape(cell2mat(data_2_use(7, 2:2:end)),5,[]),1); % Averaging CS-Trials

data.behavior_MedPc_ITI_epochs = zeros(1,(round(length(data_2_use)/5))/2 +1); 
data.behavior_MedPc_ITI_epochs(1,1) = cell2mat(data_2_use(7, 1)); % baseline
data.behavior_MedPc_ITI_epochs(1,3:2:20) = mean(reshape(cell2mat(data_2_use(7, 3:2:end)),5,[]),1); % Averaging ITI-Trials
% 
% clear ("data_2_use")

% % Retrival
% % Choose data:
% data_2_use = data.behavior_MedPc(:,1:end);
% 
% data.behavior_MedPc_CS_epochs = zeros(1,((round(length(data_2_use))+1)/2)); 
% data.behavior_MedPc_CS_epochs(1,1) = cell2mat(data_2_use(7, 1)); % baseline
% data.behavior_MedPc_CS_epochs(1,2:end) = cell2mat(data_2_use(7, 2:2:end)); % CS-Trials
% 
% clear ("data_2_use")


%% Select data to plot

% Movement
data_2_plot_1     = data.behavior_MedPc{1,1};

% Time vector
behav_MedPc_time_v = linspace(1,length(data.behavior_MedPc{1,1}),length(data.behavior_MedPc{1,1}));

% CS indexes
cs_trial       = round(data.events{3, 1} .* parameters.sampling_rate_MedPC); 

% Freezing indexes
freezing_start = data.behavior_MedPc{2, 1}(1,:);
freezing_end   = data.behavior_MedPc{2, 1}(1,:) + (data.behavior_MedPc{2, 1}(2,:) -1);

% Freezing Percentage Exposure
% data_2_plot_2  = data.behavior_MedPc{4, 1};

% Freezing Percentage Extinction, Retrieval
%data_2_plot_2  = [data.behavior_MedPc_CS_epochs];
data_2_plot_2  = cell2mat(data.behavior_MedPc(7,:));

figure
set(gcf,'color','white')
set(gcf, 'Position', get(0, 'Screensize'));

subplot(1,2,1)
hold all
%sgtitle('Habituation')
%sgtitle('Exposure')
%sgtitle('Extinction')
sgtitle('Retrieval')
%sgtitle('Renewal')


%plot(behav_bins_time_v,data_2_plot,'Color',[0.3, 0.3, 0.3, 0.3]) % raw data
plot(behav_MedPc_time_v, movmean(data_2_plot_1,5),'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6])
plot(behav_MedPc_time_v,ones(1,length(behav_MedPc_time_v)).*5,'k--')
plot([behav_MedPc_time_v(freezing_start);behav_MedPc_time_v(freezing_end)], [ones(1,length(freezing_start)).*50;ones(1,length(freezing_end)).*50],'k-','linew', 2,'Color',[.6, 0, 0])
plot([behav_MedPc_time_v(cs_trial(:,1));behav_MedPc_time_v(cs_trial(:,2))], [ones(1,length(cs_trial(:,1))).*95;ones(1,length(cs_trial(:,2))).*95],'k-','linew', 20,'Color',[1, .4, .4])


plot(behav_MedPc_time_v(100:end-100),ones(1,length(behav_MedPc_time_v(100:end-100))).*95,'k-','linew', .5,'Color',[1, .4, .4])


xlabel('Time (s)','FontSize',12), ylabel('Movement (%)','FontSize',12)
xlim([behav_MedPc_time_v(1)-10 behav_MedPc_time_v(end)])
ylim([0 100])

legend('Movement','Threshold','Freezing','','CS-Trials','NumColumns',4,'Location','southoutside')
legend('boxoff')

subplot(1,2,2)

% Exposure
% b1 = bar(data_2_plot_2);
% b1(1).FaceColor = 'w';
% ylabel('Freezing (%)','FontSize',12)
% ylim([0 105])

% Extinction, Retrieval
plot(data_2_plot_2,'-o','linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6],'MarkerEdgeColor',[.6, 0, 0])
ylabel('Freezing (%)','FontSize',12)
xticks([1:length(data_2_plot_2)])
%xticklabels({'Baseline','CS 1','CS 2','CS 3','CS 4','CS 5'})                                     % CS-Trials
%xticklabels({'Baseline','CS 1','ITI','CS 2','ITI','CS 3','ITI','CS 4','ITI','CS 5','end'})       % CS & ITI-Trials
xticklabels({'Baseline','Blk 1','Blk 2','Blk 3','Blk 4','Blk 5','Blk 6','Blk 7','Blk 8','Blk 9'}) % Average block
xtickangle(90)
xlim([0.2 length(data_2_plot_2)+1])
ylim([-5 105])


% Clear

clear('behav_bins_time_v','cs_trial','data_2_plot_1','data_2_plot_2','freezing_end','freezing_start')

%% Save

ms = 1;
%newStr = regexprep(files.id.name,'.mat','_');
newStr1 = files.id(ms).name(1:end-16);
path = files.FilesLoaded{1, 1}.folder;
%path = '/Users/flavio/Desktop';

name = strcat(path,'/',newStr1,'_behavior_extinction');

saveas(gcf,name,'png')

close all

clear('name','newStr1','path')


%% last update 18/01/2024 - 18:41
%  listening: Mogwai - Helicon 1

