function flys=Find_Ang_Vel_V2(flys)
Fs=250;
for i=1:length(flys)
    Ang_velUF=diff(flys(i).Motion_NoSaccade_Zeroed)/(1/Fs);
    [b, a] = butter(5, 7/(Fs/2),'low');
    Fil_AV = filtfilt(b, a, Ang_velUF);
    flys(i).Ang_veluf=Ang_velUF;
    flys(i).Ang_vel=Fil_AV;
    flys(i).Delta_Ang_Vel=(flys(i).Motion_NoSaccade_Zeroed(end)-flys(i).Motion_NoSaccade_Zeroed(1))/(length(flys(i).Motion_NoSaccade_Zeroed)/Fs);
    flys(i).ang_vel_mean=mean(flys(i).Ang_vel);
    flys(i).ang_vel_med=median(flys(i).Ang_vel);
end
figure
plot(flys(1).Ang_vel)
hold on
plot(flys(1).Ang_veluf)
for i=1:length(flys)
    ang_vel_test(i)=flys(i).Delta_Ang_Vel;
    ang_vel_mean(i)=flys(i).ang_vel_mean;
    ang_vel_med(i)=flys(i).ang_vel_med;
    avg_ang_disp(i)=flys(i).Motion_NoSaccade_Zeroed(end);
end
[a b]=max(abs(avg_ang_disp))
figure
boxplot(abs(ang_vel_test),ones(1,length(ang_vel_test)))
hold on
boxplot(abs(ang_vel_mean),2*ones(1,length(ang_vel_mean)),'Colors','c')
hold on
boxplot(abs(ang_vel_med),3*ones(1,length(ang_vel_med)),'Colors','g')
figure
boxplot(abs(avg_ang_disp))
title('Angular Displacement')
figure
boxplot([abs(ang_vel_test) abs(ang_vel_mean) abs(ang_vel_med)]  ,[ones(1,length(ang_vel_test)) 2*ones(1,length(ang_vel_mean)) 3*ones(1,length(ang_vel_med))])
title('Box plots for velocity by disp, mean, and median')

[h,p_del_mean]=ttest2(ang_vel_test,ang_vel_mean,'Vartype','unequal','Tail','both');
[h,p_med_mean]=ttest2(ang_vel_med,ang_vel_mean,'Vartype','unequal','Tail','both');
[h,p_del_med]=ttest2(ang_vel_test,ang_vel_med,'Vartype','unequal','Tail','both');
figure 
plot([p_del_mean p_med_mean p_del_med],'*')
title('p-value for the different methods to study the angular velocity')
ylim([-0.1 1.2])

end



