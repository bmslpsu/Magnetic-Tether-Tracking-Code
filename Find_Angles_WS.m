function [angles]=Find_Angles_WS(files, root,dirVid,dirXY,dirThresh,dirAng)
global patstype

a=length(string(files));
if a==0
    disp('No data can be loaded')
end

for j=1:length(string(files))
    v = VideoReader([dirVid files{j}]); %loads the video file
    %%
    i=1;
    while hasFrame(v)
        video_data(:,:,i)=readFrame(v);
        i = i + 1;
    end
    clear v
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
    clear A3 A4 A5
    %% finds the contour for the first proccessed image
    [C,h]=imcontour(Video_Modified(:,:,1),1);
    %% find the distance to centroid for each frame
    for kk=1:i-1
        [C]=contourc(double(Video_Modified(:,:,kk)), [1 1]);
        del_pos=sqrt((C(1,2:end)-Center_User(1)).^2+ (C(2,2:end)-Center_User(2)).^2);
        [dist, location_max]=max(del_pos);
        Max_XY(kk,1)=C(1,location_max+1);
        Max_XY(kk,2)=C(2,location_max+1);
        angles(kk)=atan2(Max_XY(kk,2)-Center_User(2),Max_XY(kk,1)-Center_User(1)); %angle in rads
        figure(1); plot(Max_XY(kk,1),Max_XY(kk,2),'*')
        hold on
        plot(C(1,2:end),C(2,2:end)); plot(Center_User(1),Center_User(2),'*')
        xlim([-100 300])
        ylim([-100 300])
        hold off
        
    end
    temp = strsplit(char(files(j)), {'_','.'});
    Angle_unw = unwrap(angles);
    Fly_Struct(j).PatternTypes=patstype;
    Fly_Struct(j).FlyNumber=temp{2};
    Fly_Struct(j).TrialNumber=temp{4};
    Fly_Struct(j).Unf_Angles=Angle_unw;
    Fly_Struct(j).COR=Center_User;
    clear Video_Cropped
end
end


