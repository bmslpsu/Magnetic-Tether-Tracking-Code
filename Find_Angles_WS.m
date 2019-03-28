function [Fly_Struct]=Find_Angles_WS(files, root,dirVid,dirXY,dirThresh,dirAng, debug)
global patstype

a=length(string(files));
if a==0
    disp('No data can be loaded')
end

for j=1:length(string(files))
    disp(['reading video ' num2str(j)]);
    if a>1
        v = VideoReader([dirVid files{j}]); %loads the video file
    elseif a==1
        v = VideoReader([dirVid files]);
    else
    end
    i=1;
    % gets the frames of the video
    while hasFrame(v)
        video_data(:,:,i)=readFrame(v);
        i = i + 1;
        if i==1000 || i==2000 || i==4000
            disp('Video is still being loaded')
        end
    end
    disp('Done reading video')
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
    disp('Finished video cropping')
    clear video_data
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
     clear A3
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
        
        [dist, index]=maxk(del_pos,15);
     
        dist_x= C(1,index+1);
        dist_y=C(2,index+1);
        
        Max_XY(1)=mean(dist_x);
        Max_XY(2)=mean(dist_y);
        angles(kk)=atan2(Max_XY(2)-Center_User(2),Max_XY(1)-Center_User(1)); %angle in rads
        
        % statement used to detect if the algorithm has located the head
        % instead of the body
        if kk>1
            if abs(angles(kk)-angles(kk-1))>pi/3
                angles(kk)=angles(kk)-pi;
            end
        end
        
        if debug==1 %allows us to see the point being tracked within the image
            figure(1); plot(Max_XY(1),Max_XY(2),'*')
            hold on
            plot(C(1,2:end),C(2,2:end)); plot(Center_User(1),Center_User(2),'*')
            xlim([-100 300])
            ylim([-100 300])
            hold off
        end
        
    end
    if a>1
        temp = strsplit(char(files(j)), {'_','.'});
    elseif a==1
        temp= strsplit(char(files), {'_','.'});
    end
    Angle_unw = unwrap(angles*180/pi);
    Fly_Struct(j).PatternTypes=patstype;
    Fly_Struct(j).FlyNumber=temp{2};
    Fly_Struct(j).TrialNumber=temp{4};
    Fly_Struct(j).Unf_Angles=Angle_unw;
    Fly_Struct(j).COR=Center_User;
    clear Video_Cropped
    save([dirAng ['\Fly_' temp{2} 'trial_' temp{4} '.mat']],'Angle_unw')
    
end
end


