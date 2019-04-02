%% code is used to place all trials in a structure for easier analysis
clc
clear all
%% initialize save locations
root = 'J:\Data for large Arena\';
dirXY = [root 'X_Y3.75\'];
dirCent = [root 'Centroid3.75\'];
dirThresh = [root 'threshold3.75\'];
dirVid = [root 'vid3.75\'];
dirAng=[root 'angles3.75\'];
patstype='3.75';
%% load the file names
[files, dirpath] = uigetfile({'*.mat', 'MAT-files'},...
    'Pick files for computing threshold', dirAng,'MultiSelect','on');
%% checks the length of the file
a=length(string(files));
if a==0
    disp('No data can be loaded')
end

%% save the data in a struct
for j=1:length(string(files))
    if a>1
        load([dirAng files{j}]);
    elseif a==1
        load([dirAng files]);
    else
    end
    
    if a>1
        temp = strsplit(char(files(j)), {'_','.'});
    elseif a==1
        temp= strsplit(char(files), {'_','.'});
    end
    Fly_Struct(j).PatternTypes=patstype;
    Fly_Struct(j).FlyNumber=temp{2};
    Fly_Struct(j).TrialNumber=temp{4};
    Fly_Struct(j).Unf_Angles=Angle_unw;
end