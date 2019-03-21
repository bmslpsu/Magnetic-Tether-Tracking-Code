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
%% Image Work modifies the video for data extraction
se= strel('square',7);
se2= strel('square',8);
A3=imdilate(video_data,se,'full');
A3=imerode(A3,se);
A4=imbinarize(A3,0.25);
A5=uint8(255*imcomplement(A4));
A5=imdilate(A5,se,'full');

%% allows the user to do the proper cropping. This decreases file size used for analysis
figure(1);
[Cropped_First_Frame, rect]=imcrop(First_Frame); %lets the user select the ROI
title('Select the region of interest for cropping')
close (1)
Video_Cropped=imcrop(A5,rect); %crops the entire video so they have the same size

%% finds the contour for the first proccessed image
[C,h]=imcontour(A5(:,:,1),1);

%% finds the centroid of the contour
centr=[mean(C(1,2:end)) mean(C(2,2:end))];
hold on
plot(centr(1),centr(2),'*')

%% find the distance to centroid 
del_pos=sqrt((C(1,2:end)-centr(1)).^2+ (C(2,2:end)-centr(2)).^2);
[dist, location_max]=max(del_pos);
ref_XY=[C(1,location_max+1),C(2,location_max+1)];
plot(ref_XY(1),ref_XY(2),'*')
%% filtering: might not be needed
% y = sgolayfilt(C(2,2:end),2,5);
% plot(C(1,2:end),y)



%% get the center of rotation from the user
figure(1) ; clf ; imshow(Cropped_First_Frame); title('Click on the center of rotation')
Center_User=ginput(1);
close (1)
angle=atan2(ref_XY(2)-Center_User(2),ref_XY(1)-Center_User(1)); %angle in rads

