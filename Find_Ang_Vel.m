function flys=Find_Ang_Vel(flys)
Fs=60;
for i=1:length(flys)
    Ang_velUF=diff(flys(i).data)/(1/Fs);
    [b, a] = butter(5, 7/(Fs/2),'low');
    Fil_AV = filtfilt(b, a, Ang_velUF);
    flys(i).Ang_veluf=Ang_velUF;
    flys(i).Ang_vel=Fil_AV;
    flys(i).Delta_Ang_Vel=(flys(i).data(end)-flys(i).data(1))/(length(flys(i).data)/60);
    flys(i).ang_vel_mean=mean(flys(i).Ang_vel);
    flys(i).ang_vel_med=median(flys(i).Ang_vel);
end
figure
plot(flys(1).Ang_vel)
hold on
plot(flys(1).Ang_veluf)
ang_vel0=[];ang_vel1=[];ang_vel2=[];ang_vel3=[];ang_vel4=[];ang_vel5=[];ang_vel6=[];ang_vel7=[];
ang_vel0_mean=[];ang_vel1_mean=[];ang_vel2_mean=[];ang_vel3_mean=[];ang_vel4_mean=[];
ang_vel5_mean=[];ang_vel6_mean=[];ang_vel7_mean=[];
for i=1:length(flys)
    if flys(i).position==0
        ang_vel0=[ang_vel0 flys(i).Delta_Ang_Vel];
        ang_vel0_mean=[ang_vel0_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==1
        ang_vel1=[ang_vel1 flys(i).Delta_Ang_Vel];
        ang_vel1_mean=[ang_vel1_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==2
        ang_vel2=[ang_vel2 flys(i).Delta_Ang_Vel];
        ang_vel2_mean=[ang_vel2_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==3
        ang_vel3=[ang_vel3 flys(i).Delta_Ang_Vel];
        ang_vel3_mean=[ang_vel3_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==4
        ang_vel4=[ang_vel4 flys(i).Delta_Ang_Vel];
        ang_vel4_mean=[ang_vel4_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==5
        ang_vel5=[ang_vel5 flys(i).Delta_Ang_Vel];
        ang_vel5_mean=[ang_vel5_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==6
        ang_vel6=[ang_vel6 flys(i).Delta_Ang_Vel];
        ang_vel6_mean=[ang_vel6_mean flys(i).ang_vel_mean];
        
    elseif flys(i).position==7
        ang_vel7=[ang_vel7 flys(i).Delta_Ang_Vel];
        ang_vel7_mean=[ang_vel7_mean flys(i).ang_vel_mean];
    end
    ang_vel_test(i)=flys(i).Delta_Ang_Vel;
    pos_test(i)=flys(i).position;
    ang_vel_mean(i)=flys(i).ang_vel_mean;
    ang_vel_med(i)=flys(i).ang_vel_med;
    avg_ang_disp(i)=flys(i).data(end);
    
end
length(ang_vel_mean)
figure
boxplot(abs(ang_vel_test),pos_test)
hold on
boxplot(abs(ang_vel_mean),pos_test)
hold on
boxplot(abs(ang_vel_med),pos_test,'Colors','g')
figure
boxplot(abs(avg_ang_disp),pos_test)
%218
