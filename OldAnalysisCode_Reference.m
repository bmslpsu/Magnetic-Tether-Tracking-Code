clear all
clc
close all
%%

root = 'E:\WaelsCodes - Copy\MagnoScript_PaperPatterns\';
global dirXY dirCent dirThresh dirVid dirAng patstype
num=3;
Choose_Patterns_to_Analyze(root,num)

%% loads the file names into matlab
[files, dirpath] = uigetfile({'*.mat', 'MAT-files'},...
    'Pick files for computing threshold', dirVid,'MultiSelect','on');
%% Load data and find angles
Fly_Struct =Find_Angles(dirVid, dirXY,dirAng, files, patstype,num);


%% Filter Data and remove saccades
Fs=160; %camera frame rate

for i=1:length(Fly_Struct)
    
    [AbsDwinN, Fwin]=Get_Abs_Graphs(Fs,Fly_Struct(i));
    motion= Remove_Saccads(AbsDwinN,Fwin,0); %find the rotation of the fly after filtering out saccades
    Fly_Struct(i).Motion_NoSaccade=motion;
    Fly_Struct(i).Fil_Angles=Fwin;
end

%% modify data
% this part of the code makes it all the same length for easy plotting and
% analysis
[Fly_Struct, min_length]=ZeroAndCropData(Fly_Struct);
time=(1:min_length)/160;
save('Fly_struct_noSaccades_3.75.mat','Fly_Struct')
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
title('Angular Displacement for each individual fly with a 7.5 degree pattern')
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
boxplot(DISP,g,'Labels',{'7.5 cw rotation','7.5 ccw roation'})
title('Displacement in CW and ccw direction for 7.5 degree patterns')

disp_total= [disp75cw; abs(disp75ccw)];
figure
boxplot(disp_total)
title('Absolute Displacement for 7.5 degree patterns')
%% Saccade Analysis
[Fly_Struct] =Saccad_Indicies(Fly_Struct);
%% Saccade Dynamics
[Fly_Struct]=Saccad_Dynamics(Fly_Struct,time);

%% Functions
function[]=  Choose_Patterns_to_Analyze(root, num)
global dirXY dirCent dirThresh dirVid dirAng patstype
if num==1
    dirXY = [root 'X_Y\'];
    dirCent = [root 'Centroid\'];
    dirThresh = [root 'threshold\'];
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

function [Fly_Struct]=Saccad_Dynamics(Fly_Struct,time)

%this function determines the saccad dynamics. the duration and average
%anglar velocity during a saccad

for i=1:length(Fly_Struct)
    Sac_start=Fly_Struct(i).saccad_Start;
    Sac_end=Fly_Struct(i).saccad_End;
    Ang_vel=diff(Fly_Struct(i).Fil_Angles)*160; %angular velocity of the saccad
    %all variables need to be reintialized to empty. This is to avoid any
    %error that will occur if the next trial has less saccads than a
    %previous trial.
    deltaT=zeros(length(Sac_start),1);
    deltaPos=zeros(length(Sac_start),1);
    Saccad_vel=zeros(length(Sac_start),1);
    Saccad_peak_vel=zeros(length(Sac_start),1);
    
    for j=1:length(Sac_end)
        deltaT(j)=time(Sac_end(j))-time(Sac_start(j)); %Duration of each saccad
        deltaPos(j)=Fly_Struct(i).Fil_Angles(Sac_end(j))-Fly_Struct(i).Fil_Angles(Sac_start(j)); %displacement of each saccad
        
        a_vel=[]; %angular velocity for each point in the saccad
        for w=Sac_start:Sac_end
            a_vel=[a_vel Ang_vel(w)];
        end
        Saccad_vel(j)=mean(a_vel);
        Saccad_peak_vel(j)=max(abs(a_vel));
    end
    
    Fly_Struct(i).count_saccad=size(Sac_end);
    Fly_Struct(i).disp_saccad=deltaPos;
    Fly_Struct(i).duration_saccad=deltaT;
    Fly_Struct(i).velocity_saccad=Saccad_vel;
    Fly_Struct(i).peak_velocity_saccad=Saccad_peak_vel;
    
    
end %for loop 1

end

function [Fly_Struct] = Saccad_Indicies(Fly_Struct)
[Fly_Struct] =Find_Saccad_Location(Fly_Struct); %debug note: this function works. but tolerence needs to be modified for different frame rates and setups
%for jeans code Find_Saccad_Location is actually used for the purpose of
%analysis. The rest of code is used as portion of another analysis
%performed by me
end


function [flys] =Find_Saccad_Location(flys)
%this function filters out the saccads using the derivative of the angular
for i =1:length(flys)
    
    Fwin=flys(i).Fil_Angles; %Filtered Angles with saccades
    Dwin = diff(Fwin)/(1/160); %diff takes an array and gives approximate derivative
    AbsDwin = abs(Dwin);
    Mvel = nanmedian(AbsDwin./0.6745); %this is done by trail and fly number, perhaphs it can be done with file name
    STDvel = nanstd(AbsDwin);
    T = Mvel + 2.5*STDvel; %tolerence at which movement is filtered out
    AbsDwinN = AbsDwin;%angualr velocity deg/sec
    AbsDwinN(AbsDwinN < T) = 0; %filteres other movements and keeps saccads
    flys(i).AbsDerWinN=AbsDwinN;
    flys(i).AbsDerWin=AbsDwin;
    % returns the index of the peak vel of the saccad according to jeans code
    % find local maxima
    S=[];
    Es=[];
    [pks, pind] = findpeaks(AbsDwinN,'MINPEAKDISTANCE',3);
    if max(pks)<45 %eliminates high speed movements which arent saccades
        
    else
        % if the max velocity of saccads is above the saccad threshold
        keep_indicies=find(pks>45); %removes the indicies in which the saccad velocity is too low
        pks=pks(keep_indicies);
        pind=pind(keep_indicies);
        pindN = pind;
        pksN=pks;
        for jj=1:length(pindN)
            %this is used to find the start and end index of the saccad
            if pindN(jj)>=length(Fwin)
            else
                if pindN(jj)<=100 ||  pindN(jj)>length(Fwin)-100
                    inc=0;
                else
                    inc=100;
                end
                saccade_start=pindN(jj)-(inc-find(AbsDwin(pindN(jj)-inc:pindN(jj)) <= pksN(jj)/4,1,'last'))-1;
                saccade_end=pindN(jj)+(find(AbsDwin(pindN(jj):end) <= pksN(jj)/4,1,'first'))-1;
                
                if isempty(saccade_start) || isempty(saccade_end)
                    saccade_start=0;
                    saccade_end=0;
                else
                    S =[S saccade_start]; %this returns the index of data left of the peak
                    Es =[Es saccade_end]; %find index of first point to the right of the peak that meets requirements of 1/4 peak velocity
                end
            end
        end
        plot(AbsDwinN)
        hold on
    end
    flys(i).saccad_Start=S;
    flys(i).saccad_End=Es;
end

end %function end

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

function []= plot_direction(cw, ccw,time)
%% subplot cw and ccw motion together

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

function [AbsDwinN, Fwin]=Get_Abs_Graphs(Fs,Fly_Struct)
% remove jumps in angle
wun = Fly_Struct.Unf_Angles;

% filter data
[b, a] = butter(3, 3.5/(Fs/2),'low');
Fwin = filtfilt(b, a, wun);


% get angular velocity in deg/s
Dwin = diff(Fwin)/(1/Fs); %diff takes an array and gives approximate derivative
AbsDwin = abs(Dwin(1:end-15)); % take abs. value and remove spike at end

%AbsDwin(AbsDwin > 3000) = NaN; % ignore peak vel > 3000

Mvel = nanmedian(AbsDwin./0.6745); %this is done by trail and fly number, perhaphs it can be done with file name
STDvel = nanstd(AbsDwin);

% Mvel2=nanmedian(AbsDwin./0.6745); %here data is sorted based on file names as an array
% STDvel= nanstd(AbsDwin);

T = Mvel + 2.5*STDvel; %tolerence at which movement is filtered out

AbsDwinN = AbsDwin;%angualr velocity deg/sec
AbsDwinN(AbsDwinN < T) = 0; %filteres other movements and keeps saccads
end

function Fly_Struct =Find_Angles(dirVid, dirXY,dirAng, files, patstype,num)

% Load the data and analyze the angles
i=1;
a=length(string(files));
while i<=a
    if a>1
        fname = [dirVid char(files(i))];
        temp = strsplit(char(files(i)), {'_','.'});
    elseif a==1
        fname = [dirVid char(files)];
        temp = strsplit(char(files), {'_','.'});
    else
        disp('No data can be loaded')
    end
    load(fname)
    [Body_Angles,centerPoint]=HeadTracker(vid_data,time_VWN,1,50,1);
    
    Angle_unw = unwrap(Body_Angles);
    Fly_Struct(i).PatternTypes=patstype;
    Fly_Struct(i).FlyNumber=temp{2};
    Fly_Struct(i).TrialNumber=temp{4};
    Fly_Struct(i).Unf_Angles=Angle_unw;
    Fly_Struct(i).time=time_VWN;
    Fly_Struct(i).COR=centerPoint;
    %loop used in order to change ROI and center of rotation
    while true
        R = input('Is the analysis satisfactory? [y/n]: ', 's');
        try
            R = validatestring( R, { 'y', 'n' } );
            switch R
                case 'y'
                    disp('saving angles, time vector, and center of rotation');
                    if a>1
                        save([dirAng files{i} '.mat'], 'Body_Angles', 'time_VWN');
                        save([dirXY files{i} '.mat'], 'centerPoint');
                        disp('finished save')
                    elseif a==1
                        save([dirAng files '.mat'], 'Body_Angles', 'time_VWN');
                        save([dirXY files '.mat'], 'centerPoint');
                        disp('finished save')
                    else
                        disp('File cannot be saved')
                    end
                    i=i+1;
                    
                case 'n'
                    disp('Trial will be repeated')
            end
        catch
            warning('Did not understand input. Try again')
            continue
        end
        break
    end
    clear R
end
if num==1
    save('flys_Data_7.5only.mat','Fly_Struct') %saves the data in a single structure
elseif num==2
    save('flys_Data_7.5.mat','Fly_Struct') %saves the data in a single structure
elseif num==3
    save('flys_Data_3.75.mat','Fly_Struct') %saves the data in a single structure
elseif num==3
    save('flys_Data_2.5.mat','Fly_Struct') %saves the data in a single structure
end

end

