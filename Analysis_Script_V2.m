clear all
clc
close all
%% notes on errors/bugs in the code
% the value of the angles must be cleared when moving to another fly.
% During a test run i noticed many trials had similar motion towards the
% end. While the data was cropped either ways, the function should still be
% modified for this error.
%%

root = 'J:\Data for large Arena\';
global dirXY dirCent dirThresh dirVid dirAng patstype
num=3; %2 for 7.5 3 for 3.75
Choose_Patterns_to_Analyze(root,num)

%% loads the file names into matlab
[files, dirpath] = uigetfile({'*.avi', 'MAT-files'},...
    'Pick files for computing threshold', dirVid,'MultiSelect','on');


%% find the yaw motion of each fly
Fly_Struct=Find_Angles_WS(files, root,dirVid,dirXY,dirThresh,dirAng,1);

%% loads the required data into the system
% this for analysis of data code not video, so you wont need the first part
% if you lareayd have the angles
if num==2
    load('flies_7.5.mat');
    Fly_Struct=Fly_Struct2;
elseif num==3
    load('flies_3.75.mat');
end
%% Filter Data and remove saccades
Fs=250; %camera frame rate

for i=1:length(Fly_Struct)
    
    [AbsDwinN, Fwin]=Get_Abs_Graphs(Fs,Fly_Struct(i));
    motion= Remove_Saccads(AbsDwinN,Fwin,0); %find the rotation of the fly after filtering out saccades
    Fly_Struct(i).Motion_NoSaccade=motion;
    Fly_Struct(i).Fil_Angles=Fwin;
    plot(motion)
    title('Motion of each fly without saccades after filtering')
    hold on
end

%% modify data
% this part of the code makes it all the same length for easy plotting and
% analysis
[Fly_Struct, min_length]=ZeroAndCropData(Fly_Struct);
time=(1:min_length)/250;
figure
for i=1:length(Fly_Struct)
    plot(time, Fly_Struct(i).Motion_NoSaccade_Zeroed)
    hold on
end

%% Find the direction count for each trial
% data in count is formated in the following order :count =[cw ccw no_rotatation]
[cw, ccw, no_rot, count]=Find_Directionv3(Fly_Struct);
bar(count)
%%
plot_direction(cw, ccw,time);
%% Means and averages for CW and CCW

figure
[medcw, gmedcw, cwflys_list]= Get_mean_avg(cw,min_length);
plot(time,gmedcw)
ylim([-100 200])
title('cw med and avg for pos')

figure
[medccw, gmedccw,  ccwflys_list]= Get_mean_avg(ccw,min_length);
plot(time,gmedccw)
ylim([-200 100])
title('ccw med and avg for pos')
%% plots position on its own
figure
for i=1:length(Fly_Struct)
    plot(time,Fly_Struct(i).Motion_NoSaccade_Zeroed)
    hold on
end
if num ==2
    title('Angular Displacement for each individual fly with a 7.5 degree pattern')
elseif num==3
    title('Angular Displacement for each individual fly with a 3.75 degree pattern')
end
xlabel('time')
ylabel('Displacement Degrees')


%% box plot for displacement
clear disp75cw disp75ccw
for j=1:length(cw)
    disp75cw(j)=cw(j).Motion_NoSaccade_Zeroed(end);
end
for j=1:length(ccw)
    disp75ccw(j)=ccw(j).Motion_NoSaccade_Zeroed(end);
end
disp75cw=transpose(disp75cw);
disp75ccw=transpose(disp75ccw);
g = [zeros(length(disp75cw), 1); ones(length(disp75ccw), 1)];
DISP=[disp75cw; disp75ccw];
figure
if num==2
    boxplot(DISP,g,'Labels',{'7.5 cw rotation','7.5 ccw roation'})
    title('Displacement in CW and ccw direction for 7.5 degree patterns')
elseif num==3
    boxplot(DISP,g,'Labels',{'3.75 cw rotation','3.75 ccw roation'})
    title('Displacement in CW and ccw direction for 3.75 degree patterns')
end

disp_total= [disp75cw; abs(disp75ccw)];
figure
boxplot(disp_total)
if num==2
    title('Absolute Displacement for 7.5 degree patterns')
elseif num==3
    title('Absolute Displacement for 3.75 degree patterns')
end

%% Functions
%---------------------------------------------------------------------
function[]=  Choose_Patterns_to_Analyze(root, num)
global dirXY dirCent dirThresh dirVid dirAng patstype
if num==1
    dirXY = [root 'X_Y\'];
    dirVid = [root 'vid\'];
    dirAng=[root 'angles\'];
    patstype='7.5';
elseif num==2
    dirXY = [root 'X_Y7.5\'];
    dirCent = [root 'Centroid7.5\'];
    dirThresh = [root 'threshold7.5\'];
    dirVid = [root 'vid7.5\'];
    dirAng=[root 'angles7.5\'];
    patstype='7.5';
elseif num==3
    dirXY = [root 'X_Y3.75\'];
    dirCent = [root 'Centroid3.75\'];
    dirThresh = [root 'threshold3.75\'];
    dirVid = [root 'vid3.75\'];
    dirAng=[root 'angles3.75\'];
    patstype='3.75';
elseif num==4
    dirXY = [root 'X_Y2.5\'];
    dirCent = [root 'Centroid2.5\'];
    dirThresh = [root 'threshold2.5\'];
    dirVid = [root 'vid2.5\'];
    dirAng=[root 'angles2.5\'];
    patstype='2.5';
else
    disp('Input value for trial selection is out of range')
    
end
end

function [AbsDwinN, Fwin]=Get_Abs_Graphs(Fs,Fly_Struct)
% remove jumps in angle
wun = Fly_Struct.Unf_Angles;

% filter data
[b, a] = butter(3, 3.5/(Fs/2),'low');
Fwin1 = filtfilt(b, a, wun);
[b1, a1] = butter(3, 3.5/50,'low');
Fwin = filtfilt(b, a, Fwin1);


% get angular velocity in deg/s
Dwin = diff(Fwin)/(1/Fs); %diff takes an array and gives approximate derivative
AbsDwin = abs(Dwin(1:end-15)); % take abs. value and remove spike at end

%AbsDwin(AbsDwin > 3000) = NaN; % ignore peak vel > 3000

Mvel = nanmedian(AbsDwin./0.6745); %this is done by trail and fly number, perhaphs it can be done with file name
STDvel = nanstd(AbsDwin);

% Mvel2=nanmedian(AbsDwin./0.6745); %here data is sorted based on file names as an array
% STDvel= nanstd(AbsDwin);

T = Mvel + 1.5*STDvel; %tolerence at which movement is filtered out

AbsDwinN = AbsDwin;%angualr velocity deg/sec
AbsDwinN(AbsDwinN < T) = 0; %filteres other movements and keeps saccads
end

function [motion] = Remove_Saccads(AbsDwinN,Fwin,debug)

% the debug option allows to compare the plots with saccads and no sacccads for better tuning the function
% the value of the tolerence still has to be tuned to improve results
%this portion gives saccad same value
maxDwin=max(AbsDwinN);
FilteredDwin=transpose(zeros(length(AbsDwinN),1));
for j=1:length(AbsDwinN)
    if AbsDwinN(j)>0
        FilteredDwin(j)=maxDwin;
    end
end
%--------------------------------------------------------------------
%--------------------------------------------------------------------
%--------------------------------------------------------------------
%---------------------------------------------------------------------


%now find the index at which the saccad starts and end
non_saccad_values=transpose(zeros(length(FilteredDwin),1));%contains index of non-saccad movement
for j=1:length(AbsDwinN)
    if FilteredDwin(j)==0
        non_saccad_values(j)=j;
    end
end

%script to remove saccads from fly's heading
motion=[];
for j=1:length(non_saccad_values)
    if non_saccad_values(j)~=0
        motion=[motion Fwin(non_saccad_values(j))];
    end
end
tol=1.1;
for j=1:length(motion)-1
    if abs(motion(j+1)-motion(j))>tol
        motion(j+1:end)=motion(j+1:end)-(motion(j+1)-motion(j));
    end
end
% this portion is used to compare the input and output data of the saccad
% filter. normally the tolerence has to be tuned to obtain optimal results
if debug==1
    plot(motion)
    hold on
    plot(Fwin)
    
end
end

function [Fly_Struct,min_length]=ZeroAndCropData(Fly_Struct)
for i=1:length(Fly_Struct)
    lengths_data(i)=length(Fly_Struct(i).Motion_NoSaccade);
    min_length=min(lengths_data); %finds array with smallest length
end
for j=1:length(Fly_Struct)
    Fly_Struct(j).Motion_NoSaccade=Fly_Struct(j).Motion_NoSaccade(1: min_length); %cuts all arrays to have same length
end
for i=1:length(Fly_Struct)
    Fly_Struct(i).Motion_NoSaccade_Zeroed=Fly_Struct(i).Motion_NoSaccade-Fly_Struct(i).Motion_NoSaccade(1); %now all data starts from 0
end
end

function []= plot_direction(cw, ccw,time)
% subplot cw and ccw motion together

figure
size(time)
for i=1:length(cw)
    plot(time,cw(i).Motion_NoSaccade_Zeroed)
    hold on
end
ylim([-30 2000])

title('Fly Displacement for Various Positions in Clockwise Direction')
figure
for i=1:length(ccw)
    plot(time,ccw(i).Motion_NoSaccade_Zeroed)
    hold on
end
ylim([-2000 200])
title('Fly Displacement for Various Positions in Counter-Clockwise Direction')
end

function [med, grand_med,fly_list]= Get_mean_avg(flys,min_length)
%fTHIS function finds the average and median for each fly and then finds
%the grand mean and average based on those
%average of averages and mean of means fts

% Finds the median for position x and grand median
fly_list=[];
for j=1:29
    AA=[];
    for i=1:length(flys)
        if str2num(flys(i).FlyNumber)==j
            AA=[AA; flys(i).Motion_NoSaccade_Zeroed];
        end
    end
    [a b]= size(AA);
    if a==1 %seems some flys have only a single trail for position 1 this eliminates the need to find a mean for these
        med(j,:)=AA;
        fly_list=[fly_list j];
    elseif isempty(AA)==1
        med(j,:)=nan*ones(1,min_length);
    else
        med(j,:)=median(AA);
        fly_list=[fly_list j];
    end
end
%filter out the nan values for the median
med(any(isnan(med),2),:)=[];
grand_med=median(med);
end

function [cw, ccw, no_rot, count]=Find_Directionv3(flys)
cw=[];
ccw=[];
no_rot=[];
for i=1:length(flys)
    avg=flys(i).Motion_NoSaccade_Zeroed(end)-flys(i).Motion_NoSaccade_Zeroed(1);
    if avg>0
        cw=[cw flys(i)];
    elseif avg<0
        ccw=[ccw flys(i)];
    else
        no_rot=[no_rot flys(i)];
    end
    count=[length(cw) length(ccw) length(no_rot)];
end
end
