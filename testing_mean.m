clear all;
load liam_pupil.mat

% data parameter setup
N_left = length(PupilLeft);
N_right = length(PupilRight);
Fs = 60;
dt = 1/Fs;
t = Timestamp.*0.001;       % timestamp in seconds

% setup problem set and baseline set data range
problem_set = [312 470 662 820 1012 1170];
baseline_set = [300 312 650 662 1000 1012];
for i=1:6
    k(i) = find(t >=problem_set(i),1);
    x(i) = find(t >=baseline_set(i),1);
end

% setup problem set time range
easy_period = (k(1):k(2))';
medium_period = (k(3):k(4))';
hard_period = (k(5):k(6))';

% setup baseline set time range
easy_base_period = (x(1):x(2))';
medium_base_period = (x(3):x(4))';
hard_base_period = (x(5):x(6))';

% clear missing data
clearPupilLeft = ClearMissingPupil(PupilLeft,N_left);
subplot(2,1,1)
plot(t, clearPupilLeft);
clearPupilRight = ClearMissingPupil(PupilRight,N_right);
subplot(2,1,2)
plot(t, clearPupilRight);

% find mean of both eyes
clearPupilMean = (clearPupilLeft+clearPupilRight)./2;

% Apply 4Hz low pass filter to remove jittering effect
Fc = 4;
wn = (2/Fs)*Fc;
b = fir1(12, wn);
clearPupilMean_fil = filter(b, 1, clearPupilMean);
figure;
hold on;
plot(t,clearPupilMean);
plot(t,clearPupilMean_fil);
title('Average Pupil Diameter of Both Eyes');
xlabel('Time /s');
ylabel('Pupil Diameter /mm');
legend('Raw','Filtered');
hold off;

% data extraction for problem
easy_raw = clearPupilMean_fil(easy_period(1):easy_period(end));
easy_t = t(easy_period(1):easy_period(end));
medium_raw = clearPupilMean_fil(medium_period(1):medium_period(end));
medium_t = t(medium_period(1):medium_period(end));
hard_raw = clearPupilMean_fil(hard_period(1):hard_period(end));
hard_t = t(hard_period(1):hard_period(end));

% data extraction for baseline
easy_base_raw = clearPupilMean_fil(easy_base_period(1):easy_base_period(end));
easy_base_t = t(easy_base_period(1):easy_base_period(end));
medium_base_raw = clearPupilMean_fil(medium_base_period(1):medium_base_period(end));
medium_base_t = t(medium_base_period(1):medium_base_period(end));
hard_base_raw = clearPupilMean_fil(hard_base_period(1):hard_base_period(end));
hard_base_t = t(hard_base_period(1):hard_base_period(end));


% plot raw problem sets
figure;
subplot(4,1,1)
plot(easy_t,easy_raw);
title('Easy Problem Set (Before Interpolation)');
xlim([312, 470]);

subplot(4,1,2)
plot(medium_t, medium_raw);
title('Medium Problem Set (Before Interpolation)');
xlim([662, 820]);

subplot(4,1,3)
plot(hard_t, hard_raw);
title('Hard Problem Set (Before Interpolation)');
ylabel('Pupil Diameter /mm');
xlabel('Time /s');
xlim([1012, 1170]);

subplot(4,1,4)
hold on
plot(easy_raw);
plot(medium_raw);
plot(hard_raw);
legend('Easy','Medium','Hard')
set(gca,'xtick',[])
hold off;

% linear interpolation
[easy_t_lin, index] = unique(easy_t);   % removed duplicate sample points
easy_raw_lin = interp1(easy_t_lin,easy_raw(index),easy_t_lin,'pchip');
[medium_t_lin, index] = unique(medium_t);   % removed duplicate sample points
medium_raw_lin = interp1(medium_t_lin,medium_raw(index),medium_t_lin,'pchip');
[hard_t_lin, index] = unique(hard_t);   % removed duplicate sample points
hard_raw_lin = interp1(hard_t_lin,hard_raw(index),hard_t_lin,'pchip');

% linear interpolation for baseline
[easy_base_t_lin, index] = unique(easy_base_t);   % removed duplicate sample points
easy_base_raw_lin = interp1(easy_base_t_lin,easy_base_raw(index),easy_base_t_lin,'pchip');
[medium_base_t_lin, index] = unique(medium_base_t);   % removed duplicate sample points
medium_base_raw_lin = interp1(medium_base_t_lin,medium_base_raw(index),medium_base_t_lin,'pchip');
[hard_base_t_lin, index] = unique(hard_base_t);   % removed duplicate sample points
hard_base_raw_lin = interp1(hard_base_t_lin,hard_base_raw(index),hard_base_t_lin,'pchip');

% plot raw data after interpolation
figure;
subplot(4,1,1)
plot(easy_t_lin,easy_raw_lin);
title('Easy Problem Set (After Interpolation)');
xlim([312, 470]);
ylim([2.5,4]);

subplot(4,1,2)
plot(medium_t_lin, medium_raw_lin);
title('Medium Problem Set (After Interpolation)');
xlim([662, 820]);

subplot(4,1,3)
plot(hard_t_lin, hard_raw_lin);
title('Hard Problem Set (After Interpolation)');
ylabel('Pupil Diameter /mm');
xlabel('Time /s');
xlim([1012, 1170]);

subplot(4,1,4)
hold on
plot(easy_raw_lin);
plot(medium_raw_lin);
plot(hard_raw_lin);
legend('Easy','Medium','Hard')
set(gca,'xtick',[])
hold off;

% plot baseline data after interpolation
figure;
subplot(4,1,1)
plot(easy_base_t_lin,easy_base_raw_lin);
title('easy baseline raw');

subplot(4,1,2)
plot(medium_base_t_lin, medium_base_raw_lin);
title('medium baseline raw');

subplot(4,1,3)
plot(hard_base_t_lin, hard_base_raw_lin);
title('hard baseline raw');

subplot(4,1,4)
hold on
plot(easy_base_raw_lin);
plot(medium_base_raw_lin);
plot(hard_base_raw_lin);
legend('Easy','Medium','Hard')
set(gca,'xtick',[])
hold off;

% extract problem in each problem set
temp_time = easy_t_lin;
for j=[1 3 5]
    temp_q = zeros(10,2);
    for i=1:10
        temp_start = problem_set(j) + 16*(i-1);
        temp_end = problem_set(j) + 12*i + 4*(i-1);
        temp_q(i,1) = find(temp_time >= temp_start,1);
        temp_q(i,2) = find(temp_time >= temp_end,1);
    end
    % set cases for extracting each question in particular problem set
    switch j
        case 1
            easy_Ts = temp_q;
            figure;
            for m=1:10
                % create cell array to store problems value
                easy_cell{m} = easy_raw_lin(easy_Ts(m,1):easy_Ts(m,2));
                hold on;
                % plot each problem on a graph
                plot(easy_raw_lin(easy_Ts(m,1):easy_Ts(m,2)));
                
            end
            legend('1','2','3','4','5','6','7','8','9','10');
            set(gca,'xtick',[]);
            ylabel('Pupil Diameter /mm');
            title('Ten Questions in Easy Problem Set');
            hold off;
            temp_time = medium_t_lin;   % move to next case
        case 3
            medium_Ts = temp_q;
            figure;
            for m=1:10
                % create cell array to store problems value
                medium_cell{m} = medium_raw_lin(medium_Ts(m,1):medium_Ts(m,2));
                hold on;
                % plot each problem on a graph
                plot(medium_raw_lin(medium_Ts(m,1):medium_Ts(m,2)));
            end
            legend('1','2','3','4','5','6','7','8','9','10');
            set(gca,'xtick',[])
            ylabel('Pupil Diameter /mm');
            title('Ten Questions in Medium Problem Set');
            hold off;
            temp_time = hard_t_lin;     % move to next case
        case 5
            hard_Ts = temp_q;
            figure;
            for m=1:10
                % create cell array to store problems value
                hard_cell{m} = hard_raw_lin(hard_Ts(m,1):hard_Ts(m,2));
                hold on;
                % plot each problem on a graph
                plot(hard_raw_lin(hard_Ts(m,1):hard_Ts(m,2)));
            end
            legend('1','2','3','4','5','6','7','8','9','10');
            set(gca,'xtick',[])
            ylabel('Pupil Diameter /mm');
            title('Ten Questions in Hard Problem Set');
            hold off;
    end
end

% average ten problems to three curve
% find average length of Ts, sample it with 60Hz to obtain Xavg
% then use interp1 to obtain sample dimension of avg version of easy_cell
% find the average of easy_cell after interp1

% find the average length of time for interpolation for easy set
for i=1:10
    N_easy_Ts(i) = floor(mean(easy_Ts(i,2) - easy_Ts(i,1) + 1));
end
N_easy_avg = floor(mean(N_easy_Ts));
easy_Ts_avg = dt*(0:N_easy_avg-1)';

% interpolate into same length of time range
for i=1:10
    N_easy_cell = length(easy_cell{i});
    temp_t = dt*(0:N_easy_cell-1)';
    easy_avgs(:,i) = interp1(temp_t,easy_cell{i},easy_Ts_avg,'pchip');
end

% find the average length of time for interpolation for medium set
for i=1:10
    N_medium_Ts(i) = floor(mean(medium_Ts(i,2) - medium_Ts(i,1) + 1));
end
N_medium_avg = floor(mean(N_medium_Ts));
medium_Ts_avg = dt*(0:N_medium_avg-1)';

% interpolate into same length of time range
for i=1:10
    N_medium_cell = length(medium_cell{i});
    temp_t = dt*(0:N_medium_cell-1)';
    medium_avgs(:,i) = interp1(temp_t,medium_cell{i},medium_Ts_avg,'pchip');
end

% find the average length of time for interpolation for hard set
for i=1:10
    N_hard_Ts(i) = floor(mean(hard_Ts(i,2) - hard_Ts(i,1) + 1));
end
N_hard_avg = floor(mean(N_hard_Ts));
hard_Ts_avg = dt*(0:N_hard_avg-1)';

% interpolate into same length of time range
for i=1:10
    N_hard_cell = length(hard_cell{i});
    temp_t = dt*(0:N_hard_cell-1)';
    hard_avgs(:,i) = interp1(temp_t,hard_cell{i},hard_Ts_avg,'pchip');
end

% calculate mean elements of each row to get mean curve of the set
easy_Avg = mean(easy_avgs,2);
medium_Avg = mean(medium_avgs,2);
hard_Avg = mean(hard_avgs,2);
figure;
hold on;
plot(easy_Ts_avg,easy_Avg);
plot(medium_Ts_avg,medium_Avg);
plot(hard_Ts_avg,hard_Avg);
legend('easy','medium','hard');
title('Averaged Problem Sets');
hold off;

% baseline subtraction
% extract 10s baseline for each problem set
% find mean of baseline for each problem set
easy_base = mean(easy_base_raw_lin);
medium_base = mean(medium_base_raw_lin);
hard_base = mean(hard_base_raw_lin);

easy_problem = easy_Avg - easy_base;
medium_problem = medium_Avg - easy_base;
hard_problem = hard_Avg - easy_base;

figure;
hold on;
plot(easy_Ts_avg,easy_problem);
plot(medium_Ts_avg,medium_problem);
plot(hard_Ts_avg,hard_problem);
legend('Easy','Medium','Hard');
title('Averaged Change in each Problem Set');
ylabel('Average Change in Pupil Diameter /mm');
xlabel('Time After Problem Onset /s')
hold off;

% find ANOVA
% extract data with 0,2,4,6,8,10,12 seconds from ten problems
% find mean of each onset time
ANOVA_set = [0 2 4 6 8 10];
for i=1:6;
    ANOVA_index(i,1) = find(easy_Ts_avg >= ANOVA_set(i),1);
    ANOVA_index(i,2) = find(medium_Ts_avg >= ANOVA_set(i),1);
    ANOVA_index(i,3) = find(hard_Ts_avg >= ANOVA_set(i),1);
end
for i=1:6
    ANOVA_data_easy(:,i) = easy_avgs(ANOVA_index(i,1),:) - easy_base;
    ANOVA_data_medium(:,i) = medium_avgs(ANOVA_index(i,2),:) - easy_base;
    ANOVA_data_hard(:,i) = hard_avgs(ANOVA_index(i,3),:) - easy_base;
end

for i=1:6
    ANOVA_median(i,1) = median(ANOVA_data_easy(:,i));
    ANOVA_median(i,2) = median(ANOVA_data_medium(:,i));
    ANOVA_median(i,3) = median(ANOVA_data_hard(:,i));
end

figure;
bar(0:2:10, ANOVA_median);
legend('Easy','Medium','Hard');
xlabel('Time from Stimuli Onset /s');
ylabel('Median of Average Change in Pupil Diameter /mm');
title('Median of Average Change in Pupil Diameter at various Time Points');

% combine ANOVA_data with other participants to find median and calcualte
% the result

% Compute ANOVA for each problem type
ANOVA_problem = (1:10)';
easy_table = table(ANOVA_problem,ANOVA_data_easy(:,1),ANOVA_data_easy(:,2),ANOVA_data_easy(:,3),ANOVA_data_easy(:,4),ANOVA_data_easy(:,5),ANOVA_data_easy(:,6),'VariableNames',{'problem','time_0','time_2','time_4','time_6','time_8','time_10'});
time_x = table([1 2 3 4 5 6]','VariableNames',{'Measurements'});
rm = fitrm(easy_table,'time_0-time_10~problem','WithinDesign',time_x);
ANOVA_result = ranova(rm)

%% export output
easy_base_liam = easy_base;
easy_avgs_liam = easy_avgs;
medium_avgs_liam = medium_avgs;
hard_avgs_liam = hard_avgs;
easy_Ts_avg_liam = easy_Ts_avg;

save('final_pupil_liam', 'easy_base_liam','easy_avgs_liam','medium_avgs_liam', 'hard_avgs_liam','easy_Ts_avg_liam');






