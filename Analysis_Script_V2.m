clear all
clc
close all
%%

root = 'E:\Data for large Arena\';
global dirXY dirCent dirThresh dirVid dirAng patstype
num=2;
Choose_Patterns_to_Analyze(root,num)

%% loads the file names into matlab
[files, dirpath] = uigetfile({'*.avi', 'MAT-files'},...
    'Pick files for computing threshold', dirVid,'MultiSelect','on');


%%
angles=Find_Angles_WS(files, root,dirVid,dirXY,dirThresh,dirAng,1);

%% Functions
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

