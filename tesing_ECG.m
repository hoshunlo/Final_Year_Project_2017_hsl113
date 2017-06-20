load liam_ECG_RR.mat

% RRI file data parameter setup
t = Raw.rr_time;
RRI = Raw.rr_val;
t_clean = Clean.rr_time;
RRI_clean  = Clean.rr_val;
RR_Fs = 4;
Fs = 1000;
N_clean = length(RRI_clean);

% interpolate 4 Hz for frequency domain
t_interp = t_clean(1)/fs:1/RR_Fs:t_clean(end);
RRI_interp = interp1(t_clean(1:end),RRI_clean(1:end),t_interp,'pchip');

% plot entire raw data
figure;
plot(t_clean,RRI_clean);
title('RR Intervals of entire ECG recordings');
ylabel('RR intervals /ms')
xlabel('Time /s');


% data parameter setup
% dt = 1/Fs;
% t = Timestamp.*0.001;       % timestamp in seconds

% setup problem set and baseline set data range
problem_set = [312 470 662 820 1012 1170];
baseline_set = [0 300 470 650 820 1000];
real_baseline = [142 300];
for i=1:6
    k(i) = find(t_interp >=problem_set(i),1);
    x(i) = find(t_interp >=baseline_set(i),1);
    if(i<=2)
        y(i) = find(t_interp >=real_baseline(i),1);
    end
end

% setup problem set time range
easy_period = (k(1):k(2))';
medium_period = (k(3):k(4))';
hard_period = (k(5):k(6))';

% setup baseline set time range
easy_base_period = (x(1):x(2))';
medium_base_period = (x(3):x(4))';
hard_base_period = (x(5):x(6))';

% setup baseline for subtraction
real_baseline_period = (y(1):y(2))';

% data extraction for problem
easy_raw = RRI_interp(easy_period(1):easy_period(end));
easy_t = t_interp(easy_period(1):easy_period(end));
medium_raw = RRI_interp(medium_period(1):medium_period(end));
medium_t = t_interp(medium_period(1):medium_period(end));
hard_raw = RRI_interp(hard_period(1):hard_period(end));
hard_t = t_interp(hard_period(1):hard_period(end));
%%
% data extraction for non interpolation
for i=1:6
    k(i) = find(t_clean >=problem_set(i),1);
    x(i) = find(t_clean >=baseline_set(i),1);
    if(i<=2)
        y(i) = find(t_clean >=real_baseline(i),1);
    end
end

% setup problem set time range
easy_period_clean = (k(1):k(2))';
medium_period_clean = (k(3):k(4))';
hard_period_clean = (k(5):k(6))';

% setup baseline set time range
easy_base_period_clean = (x(1):x(2))';
medium_base_period_clean = (x(3):x(4))';
hard_base_period_clean = (x(5):x(6))';

% setup baseline for subtraction
real_baseline_period_clean = (y(1):y(2))';

% data extraction for problem for non interp
easy_raw_clean = RRI_clean(easy_period_clean(1):easy_period_clean(end));
easy_t_clean = t_clean(easy_period_clean(1):easy_period_clean(end));
medium_raw_clean = RRI_clean(medium_period_clean(1):medium_period_clean(end));
medium_t_clean = t_clean(medium_period_clean(1):medium_period_clean(end));
hard_raw_clean = RRI_clean(hard_period_clean(1):hard_period_clean(end));
hard_t_clean = t_clean(hard_period_clean(1):hard_period_clean(end));
%%

% data extraction for baseline
easy_base_raw = RRI_interp(easy_base_period(1):easy_base_period(end));
easy_base_t = t_interp(easy_base_period(1):easy_base_period(end));
medium_base_raw = RRI_interp(medium_base_period(1):medium_base_period(end));
medium_base_t = t_interp(medium_base_period(1):medium_base_period(end));
hard_base_raw = RRI_interp(hard_base_period(1):hard_base_period(end));
hard_base_t = t_interp(hard_base_period(1):hard_base_period(end));

% data extraction for real baseline
baseline_2min = RRI_interp(real_baseline_period(1):real_baseline_period(end));
baseline_2min_t = t_interp(real_baseline_period(1):real_baseline_period(end));

% set mutual time range
problem_t = (0:length(easy_raw)-1)';
problem_t = problem_t./length(easy_raw).*158;

% set mutual time range for non interpolation
problem_t_clean = (0:length(easy_raw_clean)-1)';
problem_t_clean = problem_t_clean./length(easy_raw_clean).*158;

% plot raw problem sets
figure;
subplot(3,1,1)
plot(easy_t,easy_raw);
title('RR Intervals of Tasks in Easy Level');
xlim([312, 470]);

subplot(3,1,2)
plot(medium_t, medium_raw);
title('RR Intervals of Tasks in Medium Level');
xlim([662, 820]);

subplot(3,1,3)
plot(hard_t, hard_raw);
xlim([1012, 1170]);
title('RR Intervals of Tasks in Hard Level');
xlabel('Time /s');
ylabel('RR Intervals /ms');

figure;
hold on
plot(problem_t,easy_raw);
plot(problem_t,medium_raw);
plot(problem_t,hard_raw);
legend('Easy','Medium','Hard')
title('RR Intervals of Tasks in Different Difficulty Levels');
ylabel('RR Intervals /ms');
xlabel('Time after Problem Set Onset /s');
hold off;

% plot baseline raw data
figure;
subplot(4,1,1)
plot(easy_base_t,easy_base_raw);
title('easy raw');
xlim([0, 300]);

subplot(4,1,2)
plot(medium_base_t, medium_base_raw);
title('medium raw');
xlim([470, 650]);

subplot(4,1,3)
plot(hard_base_t, hard_base_raw);
title('hard raw');
xlim([820, 1000]);

subplot(4,1,4)
hold on
plot(easy_base_raw,'b');
plot(medium_base_raw,'r');
plot(hard_base_raw,'g');
legend('easy','medium','hard')
hold off;

% baseline subtraction
baseline_mean = mean(baseline_2min);
easy_problem = easy_raw - baseline_mean;
medium_problem = medium_raw - baseline_mean;
hard_problem = hard_raw - baseline_mean;
% easy_problem = easy_raw;
% medium_problem = medium_raw;
% hard_problem = hard_raw;
% baseline subtraction
easy_problem_clean = easy_raw_clean - baseline_mean;
medium_problem_clean = medium_raw_clean - baseline_mean;
hard_problem_clean = hard_raw_clean - baseline_mean;


figure;
subplot(3,1,1)
plot(easy_t,easy_problem);
title('easy raw after');
xlim([312, 470]);

subplot(3,1,2)
plot(medium_t, medium_problem);
title('medium raw after');
xlim([662, 820]);

subplot(3,1,3)
plot(hard_t, hard_problem);
title('hard raw after');
xlim([1012, 1170]);

figure;
hold on
plot(problem_t,easy_problem);
plot(problem_t,medium_problem);
plot(problem_t,hard_problem);
legend('easy','medium','hard')
title('RR Intervals of Tasks in Different Difficulty Levels from Baseline');
ylabel('RR Intervals from Baseline /ms');
xlabel('Time after Problem Set Onset /s');
hold off;


% compute HR
HR_easy = 60000*1/mean(easy_raw_clean)-60000*1/baseline_mean;
HR_medium = 60000*1/mean(medium_raw_clean)-60000*1/baseline_mean;
HR_hard = 60000*1/mean(hard_raw_clean)-60000*1/baseline_mean;

% compute Mean RR
RRI_Average_easy = mean(easy_problem_clean);
RRI_Average_medium = mean(medium_problem_clean);
RRI_Average_hard = mean(hard_problem_clean);

% compute sdnn
sdnn_easy = std(easy_problem_clean);
sdnn_medium = std(medium_problem_clean);
sdnn_hard = std(hard_problem_clean);

% compute sdsd
RR_diff_easy = abs(diff(easy_problem_clean));
r_easy = sqrt(sum((RR_diff_easy-mean(RR_diff_easy)).^2)/length(RR_diff_easy));
sdsd_easy = r_easy;

RR_diff_medium = abs(diff(medium_problem_clean));
r_medium = sqrt(sum((RR_diff_medium-mean(RR_diff_medium)).^2)/length(RR_diff_medium));
sdsd_medium = r_medium;

RR_diff_hard = abs(diff(hard_problem_clean));
r_hard = sqrt(sum((RR_diff_hard-mean(RR_diff_hard)).^2)/length(RR_diff_hard));
sdsd_hard = r_hard;

% compute rMSSD
rMSSD_easy = sqrt(sum(RR_diff_easy.^2)/length(RR_diff_easy));
rMSSD_medium = sqrt(sum(RR_diff_medium.^2)/length(RR_diff_medium));
rMSSD_hard = sqrt(sum(RR_diff_hard.^2)/length(RR_diff_hard));

% compute pNN50
pNN50_easy = pNNx(easy_problem_clean,50);
pNN50_medium = pNNx(medium_problem_clean,50);
pNN50_hard = pNNx(hard_problem_clean,50);

% compute pNN50
pNN50_easy = pNNx(easy_problem_clean,50);
pNN50_medium = pNNx(medium_problem_clean,50);
pNN50_hard = pNNx(hard_problem_clean,50);

% compute pNN25 
pNN25_easy = pNNx(easy_problem_clean,25);
pNN25_medium = pNNx(medium_problem_clean,25);
pNN25_hard = pNNx(hard_problem_clean,25);

%% Frequency Domain Analysis %%
% apply bandpass filter 0.04-0.4Hz
[b_bp,a_bp] = butter(4,[0.04 0.4]/(RR_Fs*0.5),'bandpass');
easy_problem_fil = filtfilt(b_bp,a_bp,easy_problem);
medium_problem_fil = filtfilt(b_bp,a_bp,medium_problem);
hard_problem_fil = filtfilt(b_bp,a_bp,hard_problem);

% setup parameters for frequency spectrum
N_problem = length(hard_problem_fil);
nfft = 2^nextpow2(N_problem);
%window = floor(RR_Fs*120); %(1/0.04)*2 = 50 approx. 1 min, so use 2 min
window = 256;
%noverlap = round(0.2*window);
noverlap = 0.5*window;

% zero padding
if(length(easy_problem_fil) < window)
    easy_problem_fil = [easy_problem_fil,zeros(1,window-length(easy_problem_fil))];
end
if(length(medium_problem_fil) < window)
    medium_problem_fil = [medium_problem_fil,zeros(1,window-length(medium_problem_fil))];
end
if(length(hard_problem_fil) < window)
    hard_problem_fil = [hard_problem_fil,zeros(1,window-length(hard_problem_fil))];
end

% setup frqeuency bands range
ULFlow = 0;
ULFhigh = 0.003;
VLFlow = 0.003;
VLFhigh = 0.04;
LFlow = 0.04;
LFhigh = 0.15;
HFlow = 0.15;
HFhigh = 0.4;

% find PSD by Welch method
% [PSD_easy,F_easy] = pwelch(easy_problem_fil,window,noverlap,nfft,RR_Fs);
% [PSD_medium,F_medium] = pwelch(medium_problem_fil,window,noverlap,nfft,RR_Fs);
% [PSD_hard,F_hard] = pwelch(hard_problem_fil,window,noverlap,nfft,RR_Fs);

% find PSD by Lomb-Scargle method 
[PSD_easy,F_easy] = plomb(easy_problem_fil, problem_t,HFhigh);
[PSD_medium,F_medium] = plomb(medium_problem_fil, problem_t,HFhigh);
[PSD_hard,F_hard] = plomb(hard_problem_fil, problem_t,HFhigh);

% find PSD by Burg method
% AR_order = 16; % trial and error on different order
% %N_problem_Burg = N_problem*(fs/HFhigh)-1;
% window_Burg = hamming(N_problem)';
% easy_problem_Burg = window_Burg.*easy_problem_fil;
% medium_problem_Burg = window_Burg.*medium_problem_fil;
% hard_problem_Burg = window_Burg.*hard_problem_fil;
% 
% [PSD_easy,F_easy] = pburg(easy_problem_fil, AR_order, N_problem, RR_Fs);
% [PSD_medium,F_medium] = pburg(medium_problem_fil,AR_order, N_problem, RR_Fs);
% [PSD_hard,F_hard] = pburg(hard_problem_fil, AR_order, N_problem, RR_Fs);

% plot PSD for all problem sets
figure;
subplot(3,1,1)
plot(F_easy,PSD_easy);
xlim([0 0.5]);
title('PSD of Easy Level by Burg Method');

subplot(3,1,2)
plot(F_medium,PSD_medium);
xlim([0 0.5]);
title('PSD of Medium Level by Burg Method');

subplot(3,1,3)
plot(F_hard,PSD_hard);
xlim([0 0.5]);
title('PSD of Hard Level by Burg Method');
xlabel('Frequency /Hz');
ylabel('Power / ms^2');

figure;
hold on;
plot(F_easy,PSD_easy);
plot(F_easy,PSD_medium);
plot(F_easy,PSD_hard);
legend('Easy','Medium','Hard');
title('PSD of Different Level by Burg Method with AR(16)');
xlabel('Frequency /Hz');
ylabel('Power / ms^2');
xlim([0 0.5]);
hold off;


% plot different AR order for Burg for hard level
% window_Burg = hamming(N_problem)';
% hard_problem_Burg = window_Burg.*hard_problem_fil;
% figure;
% for AR_order = [8 12 14 16 18]; % trial and error on different order
% 
% [PSD_hard,F_hard] = pburg(hard_problem_fil, AR_order, N_problem, RR_Fs);
% hold on;
% plot(F_hard, PSD_hard);
% end
% xlabel('Frequency /Hz');
% ylabel('Power / ms^2');
% xlim([0 0.5]);
% title('PSD of Hard Level by Burg Method with different AR order');
% legend('8th','12th','14th','16th','18th');

% find indices for VLF, LF and HF bands
VLF_index_easy = find((F_easy>VLFlow) & (F_easy<=VLFhigh));
LF_index_easy = find((F_easy>LFlow) & (F_easy<=LFhigh));
HF_index_easy = find((F_easy>HFlow) & (F_easy<=HFhigh));

VLF_index_medium = find((F_medium>VLFlow) & (F_medium<=VLFhigh));
LF_index_medium = find((F_medium>LFlow) & (F_medium<=LFhigh));
HF_index_medium = find((F_medium>HFlow) & (F_medium<=HFhigh));

VLF_index_hard = find((F_hard>VLFlow) & (F_hard<=VLFhigh));
LF_index_hard = find((F_hard>LFlow) & (F_hard<=LFhigh));
HF_index_hard = find((F_hard>HFlow) & (F_hard<=HFhigh));

% find area of particular frequency band by trapezoidal numerical
% integration
areaVLF_easy = trapz(F_easy(VLF_index_easy),PSD_easy(VLF_index_easy))*(10^6);
areaLF_easy = trapz(F_easy(LF_index_easy),PSD_easy(LF_index_easy))*(10^6);
areaHF_easy = trapz(F_easy(HF_index_easy),PSD_easy(HF_index_easy))*(10^6);

areaVLF_medium = trapz(F_medium(VLF_index_medium),PSD_medium(VLF_index_medium))*(10^6);
areaLF_medium = trapz(F_medium(LF_index_medium),PSD_medium(LF_index_medium))*(10^6);
areaHF_medium = trapz(F_medium(HF_index_medium),PSD_medium(HF_index_medium))*(10^6);

areaVLF_hard = trapz(F_hard(VLF_index_hard),PSD_hard(VLF_index_hard))*(10^6);
areaLF_hard = trapz(F_hard(LF_index_hard),PSD_hard(LF_index_hard))*(10^6);
areaHF_hard = trapz(F_hard(HF_index_hard),PSD_hard(HF_index_hard))*(10^6);

% find area of total power and LFHF band
area_totalP_easy = areaVLF_easy + areaLF_easy + areaHF_easy;
area_totalP_medium = areaVLF_medium + areaLF_medium + areaHF_medium;
area_totalP_hard = areaVLF_hard + areaLF_hard + areaHF_hard;

area_LFHF_easy = areaLF_easy + areaHF_easy;
area_LFHF_medium = areaLF_medium + areaHF_medium;
area_LFHF_hard = areaLF_hard + areaHF_hard;

% calculate the pertentage power of total power
VLF_easy = (areaVLF_easy/area_totalP_easy)*100;
LF_easy = (areaLF_easy/area_totalP_easy)*100;
HF_easy = (areaHF_easy/area_totalP_easy)*100;

VLF_medium = (areaVLF_medium/area_totalP_medium)*100;
LF_medium = (areaLF_medium/area_totalP_medium)*100;
HF_medium = (areaHF_medium/area_totalP_medium)*100;

VLF_hard = (areaVLF_hard/area_totalP_hard)*100;
LF_hard = (areaLF_hard/area_totalP_hard)*100;
HF_hard = (areaHF_hard/area_totalP_hard)*100;

% calculate normalised power of LF and HF
LF_Norm_easy = areaLF_easy/area_LFHF_easy;
HF_Norm_easy = areaHF_easy/area_LFHF_easy;

LF_Norm_medium = areaLF_medium/area_LFHF_medium;
HF_Norm_medium = areaHF_medium/area_LFHF_medium;

LF_Norm_hard = areaLF_hard/area_LFHF_hard;
HF_Norm_hard = areaHF_hard/area_LFHF_hard;

% calculate ratio of LFHF
LFHF_Ratio_easy = areaLF_easy/areaHF_easy;
LFHF_Ratio_medium = areaLF_medium/areaHF_medium;
LFHF_Ratio_hard = areaLF_hard/areaHF_hard;

%% End of Frequency Domain Method %%

%% Ouput data here
% time domain
mean_HR_easy_liam = HR_easy;
mean_HR_medium_liam = HR_medium;
mean_HR_hard_liam = HR_hard;

mean_RRI_easy_liam = RRI_Average_easy;
mean_RRI_medium_liam = RRI_Average_medium;
mean_RRI_hard_liam = RRI_Average_hard;

sdnn_easy_liam = sdnn_easy;
sdnn_medium_liam = sdnn_medium;
sdnn_hard_liam = sdnn_hard;

rMSSD_easy_liam = rMSSD_easy;
rMSSD_medium_liam = rMSSD_medium;
rMSSD_hard_liam = rMSSD_hard;

sdsd_easy_liam = sdsd_easy;
sdsd_medium_liam = sdsd_medium;
sdsd_hard_liam = sdsd_hard;

pNN50_easy_liam = pNN50_easy;
pNN50_medium_liam = pNN50_medium;
pNN50_hard_liam = pNN50_hard;

pNN25_easy_liam = pNN25_easy;
pNN25_medium_liam = pNN25_medium;
pNN25_hard_liam = pNN25_hard;

save('final_time_liam','mean_HR_easy_liam','mean_HR_medium_liam','mean_HR_hard_liam','mean_RRI_easy_liam','mean_RRI_medium_liam','mean_RRI_hard_liam','sdnn_easy_liam','sdnn_medium_liam','sdnn_hard_liam','rMSSD_easy_liam','rMSSD_medium_liam','rMSSD_hard_liam','sdsd_easy_liam','sdsd_medium_liam','sdsd_hard_liam','pNN50_easy_liam','pNN50_medium_liam','pNN50_hard_liam','pNN25_easy_liam','pNN25_medium_liam','pNN25_hard_liam');
% frequency domain
Total_Power_easy_liam = area_totalP_easy;
Total_Power_medium_liam = area_totalP_medium;
Total_Power_hard_liam = area_totalP_hard;

LF_easy_liam = areaLF_easy;
LF_medium_liam = areaLF_medium;
LF_hard_liam = areaLF_hard;

HF_easy_liam = areaHF_easy;
HF_medium_liam = areaHF_medium;
HF_hard_liam = areaHF_hard;

LF_Norm_easy_liam = LF_Norm_easy;
LF_Norm_medium_liam = LF_Norm_medium;
LF_Norm_hard_liam = LF_Norm_hard;

HF_Norm_easy_liam = HF_Norm_easy;
HF_Norm_medium_liam = HF_Norm_medium;
HF_Norm_hard_liam = HF_Norm_hard;

LFHF_Ratio_easy_liam = LFHF_Ratio_easy;
LFHF_Ratio_medium_liam = LFHF_Ratio_medium;
LFHF_Ratio_hard_liam = LFHF_Ratio_hard;

save('final_frequency_liam_lomb','Total_Power_easy_liam','Total_Power_medium_liam','Total_Power_hard_liam','LF_easy_liam','LF_medium_liam','LF_hard_liam','HF_easy_liam','HF_medium_liam','HF_hard_liam','LF_Norm_easy_liam','LF_Norm_medium_liam','LF_Norm_hard_liam','HF_Norm_easy_liam','HF_Norm_medium_liam','HF_Norm_hard_liam','LFHF_Ratio_easy_liam','LFHF_Ratio_medium_liam','LFHF_Ratio_hard_liam');

