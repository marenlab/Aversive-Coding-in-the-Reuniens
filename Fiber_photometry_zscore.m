% extract deltF/F from tables arrays (3 column) -> data organization from pMat
%  prepare data before loading from Analises_b

data_f = cell(2,size(data_onset,2));

for ii = 1:size(data_onset,2)
    data_f{1,ii} = data_onset{1,ii}(:,3)'; % onset
    data_f{2,ii} = data_offset{1,ii}(:,3)'; % offset
end

time_F = linspace(-2,2,length(data_f{1,1})); % similar for all time windows

clear ('ii','data_onset','data_offset','output')

%% Behavior all trial
%  prepare data before loading individual files

behavior_epochs = cell(2,size(behavior,2));

for jj = 1:size(behavior,2)

    temp = behavior{1, jj}.data.behavior{1,1};

    for ii = 1:size(behavior{1, 1}.data.events{1, 1},1)

        behavior_epochs{1,jj}(ii,:)  = temp(1,behavior{1, 1}.data.events{1, 1}(ii,1) - 2*30 : behavior{1, 1}.data.events{1, 1}(ii,1) + 2*30); % onset
        behavior_epochs{2,jj}(ii,:)  = temp(1,behavior{1, 1}.data.events{1, 1}(ii,2) - 2*30 : behavior{1, 1}.data.events{1, 1}(ii,2) + 2*30); % offset

    end
end


time_b = linspace(-2,2,length(behavior_epochs{1,1}));

clear('ii','jj','temp')

%% behavior mean trials

behavior_epochs_mean = cell(2,size(behavior_epochs,2));

for ii = 1:size(behavior_epochs,2)
    behavior_epochs_mean{1,ii} = mean(behavior_epochs{1,ii},1);
    behavior_epochs_mean{2,ii} = mean(behavior_epochs{2,ii},1);

end

clear('ii')

%% zscore Delta F

data_f_z = cell(size(data_f));

for jj = 1:size(data_f,1)
    
    for ii = 1:size(data_f,2)

    % delta F
    mean_pre_f{jj,ii} = mean(data_f{jj,ii}(1,1:floor(length(time_F)/2)),2,'omitnan');    
    std_pre_f{jj,ii}  = std(data_f{jj,ii}(1,1:floor(length(time_F)/2)),0,2,'omitnan');    
    data_f_z{jj,ii} = (data_f{jj,ii} - mean_pre_f{jj,ii})./std_pre_f{jj,ii};
    mean_pos_z{jj,ii} = mean(data_f_z{jj,ii}(1,floor(length(time_F)/2) + 1 : end),2,'omitnan');

    end 

    data_f_z_mean_SEM{jj,1} = mean(cat(3,data_f_z{jj,:}),3,'omitnan'); %mean
    data_f_z_mean_SEM{jj,2} = std(cat(3,data_f_z{jj,:}),[],3,'omitnan')./size(cat(3,data_f_z{jj,:}),3); %SEM
    
    data_f_mean_SEM{jj,1} = mean(cat(3,data_f{jj,:}),3,'omitnan'); %mean
    data_f_mean_SEM{jj,2} = std(cat(3,data_f{jj,:}),[],3,'omitnan')./size(cat(3,data_f_z{jj,:}),3); %SEM

end

clear ('ii','jj')

%% zscore behavior

data_behavior_z = cell(size(behavior_epochs_mean));

for jj = 1:size(behavior_epochs,1)
    
%     for ii = 1:size(behavior_epochs,2)
% 
%     mean_pre_Behavior{jj,ii} = mean(behavior_epochs_mean{jj,ii}(1,1:floor(length(time_b)/2)),2,'omitnan');    
%     std_pre_Behavior{jj,ii}  = std(behavior_epochs_mean{jj,ii}(1,1:floor(length(time_b)/2)),0,2,'omitnan');    
%     data_behavior_z{jj,ii} = (behavior_epochs_mean{jj,ii} - mean_pre_Behavior{jj,ii})./std_pre_Behavior{jj,ii};
%     mean_pos_z_behavior{jj,ii} = mean(data_behavior_z{jj,ii}(1,floor(length(time_b)/2) + 1 : floor(length(time_b)/2) + 25),2,'omitnan');
% 
%     end 

    data_behavior_z_mean_SEM{jj,1} = mean(cat(3,data_behavior_z{jj,:}),3,'omitnan'); %mean
    data_behavior_z_mean_SEM{jj,2} = std(cat(3,data_behavior_z{jj,:}),[],3,'omitnan')./size(cat(3,data_behavior_z{jj,:}),3); %SEM

    data_behavior_mean_SEM{jj,1} = mean(cat(3,behavior_epochs_mean{jj,:}),3,'omitnan'); %mean
    data_behavior_mean_SEM{jj,2} = std(cat(3,behavior_epochs_mean{jj,:}),[],3,'omitnan')./size(cat(3,behavior_epochs_mean{jj,:}),3); %SEM    
    

end

clear ('ii','jj')

%% plot

% figure
% set(gcf,'color','w');

hold on
subplot(221)
boundedline(time_F,data_f_z_mean_SEM{1, 1},data_f_z_mean_SEM{1, 2},'k') % baseline
%boundedline(time_F,data_f_z_mean_SEM{1, 1},data_f_z_mean_SEM{1, 2},'color',[.8 0 0]) % First 10 CS

xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])
ylim([-5 20])
% yticks([-8,-6,-4,-2,0,2])

xlabel('Time(s)','FontSize',18)
ylabel({'\Delta F/F'; 'Z-score'},'FontSize',18)

title('CS Onset','FontSize',18)


subplot(222)
boundedline(time_F,data_f_z_mean_SEM{2, 1},data_f_z_mean_SEM{2, 2},'k') % baseline
%boundedline(time_F,data_f_z_mean_SEM{2, 1},data_f_z_mean_SEM{2, 2},'color',[.8 0 0]) % First 10 CS

xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])
xlabel('Time(s)','FontSize',18)

ylim([-5 20])
% yticks([-2,0,2,4,6])
title('CS Offset','FontSize',18)


subplot(223)
%boundedline(time_b,data_behavior_z_mean_SEM{1,1},data_behavior_z_mean_SEM{1,2},'k') % baseline
%boundedline(time_b,data_behavior_z_mean_SEM{1,1},data_behavior_z_mean_SEM{1,2},'color',[.8 0 0]) % First 10 CS

for ii = 1:size(data_behavior_z,2)
    plot(time_b,data_behavior_z{1, ii},'color',[.6 .6 .6],'linew',1);
    %plot(time_b,data_behavior_z{1, ii},'color',[.8 0 0 .5],'linew',1);

    hold on
end

xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])

%ylim([150 1400])


xlabel('Time(s)','FontSize',18)
ylabel({'Movement' ;'Zscore'},'FontSize',18)


subplot(224)
%boundedline(time_b,data_behavior_z_mean_SEM{2,1},data_behavior_z_mean_SEM{2,2},'k') % baseline
%boundedline(time_b,data_behavior_z_mean_SEM{2,1},data_behavior_z_mean_SEM{2,2},'color',[.8 0 0]) % First 10 CS
for ii = 1:size(data_behavior_z,2)
    plot(time_b,data_behavior_z{2, ii},'color',[.6 .6 .6],'linew',1);
    %plot(time_b,data_behavior_z{2, ii},'color',[.8 0 0 .5],'linew',1);

    hold on
end

xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])

%ylim([150 1400])


xlabel('Time(s)','FontSize',18)


%% plot

figure
set(gcf,'color','w');

subplot(221)
boundedline(time_F,data_f_mean_SEM{1, 1},data_f_mean_SEM{1, 2},'k')
xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])
ylim([-.5 2.5])
% yticks([-8,-6,-4,-2,0,2])

xlabel('Time(s)','FontSize',14)
ylabel('\Delta F/F','FontSize',14)

title('Initiation','FontSize',14)



subplot(222)
boundedline(time_F,data_f_mean_SEM{2, 1},data_f_mean_SEM{2, 2},'k')
xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])
ylim([-.5 2.5])
% yticks([-2,0,2,4,6])


subplot(223)
plot(time_b,data_behavior_mean_SEM{1,1},'k')
xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])

%ylim([150 350])


xlabel('Time(s)','FontSize',14)
ylabel('Movement(a.u.)','FontSize',14)


subplot(224)
plot(time_b,data_behavior_mean_SEM{2,1},'k')
xline(0,'k--')

xlim([-2 2])
xticks([-2,-1,-0,1,2])

%ylim([150 350])


xlabel('Time(s)','FontSize',14)


