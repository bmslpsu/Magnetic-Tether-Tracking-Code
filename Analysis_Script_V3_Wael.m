clc
clear all
close all
%%
global dirXY dirCent dirThresh dirVid dirAng patstype
dirVid='E:\WaelsCodes - Copy\Large Magno 7.5 and 3.75 experiments\vid';
v = VideoReader('Fly1_Trial1_C001H001S0001.avi');
%%
i=1;
while hasFrame(v)
    video_data(:,:,i)=readFrame(v);
     i = i + 1;
end
%% take the first frame of the video
First_Frame= video_data(:,:,1);

%% thresholding
thresh=Find_Threshold(video_data);

%% allows the user to do the proper cropping. This decreases file size used for analysis
figure(1);
[Cropped_First_Frame, rect]=imcrop(First_Frame); %lets the user select the ROI
title('Select the region of interest for cropping')
close (1)
for kk=1:i-1
    Video_Cropped(:,:,kk)=imcrop(video_data(:,:,kk),rect); %crops the entire video so they have the same size
end
%% finds the centroid of the contour
figure(1) ; clf ; imshow(Cropped_First_Frame); title('Click on the center of rotation')
Center_User=ginput(1);
close (1)
%% Image Work modifies the video for data extraction
se= strel('square',10);
se2= strel('square',10);
A3=imdilate(Video_Cropped,se,'full');
A3=imerode(A3,se);
A4=imbinarize(A3,thresh);
A5=uint8(255*imcomplement(A4));
Video_Modified=imdilate(A5,se,'full');
imshow(Video_Modified(:,:,1))
%% finds the contour for the first proccessed image
[C,h]=imcontour(Video_Modified(:,:,1),1);
%% find the distance to centroid for each frame
for kk=1:i-1
    [C]=contourc(double(Video_Modified(:,:,kk)), [1 1]);
    del_pos=sqrt((C(1,2:end)-Center_User(1)).^2+ (C(2,2:end)-Center_User(2)).^2);
    [dist, location_max]=max(del_pos);
    Max_XY(kk,1)=C(1,location_max+1);
    Max_XY(kk,2)=C(2,location_max+1);
    angle(kk)=atan2(Max_XY(kk,2)-Center_User(2),Max_XY(kk,1)-Center_User(1)); %angle in rads
    figure(1); plot(Max_XY(kk,1),Max_XY(kk,2),'*')
    hold on 
    plot(C(1,2:end),C(2,2:end)); plot(Center_User(1),Center_User(2),'*')
    xlim([-100 300])
    ylim([-100 300])
    hold off

end




