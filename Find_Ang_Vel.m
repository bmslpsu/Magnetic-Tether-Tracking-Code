function Fly_Struct=Find_Ang_Vel(Fly_Struct)
Fs=160;
for i=1:length(Fly_Struct)
    Ang_velUF=diff(Fly_Struct(i).Motion_NoSaccade_Zeroed)/(1/Fs);
    [b, a] = butter(5, 15/(Fs/2),'low');
    Fil_AV = filtfilt(b, a, Ang_velUF);
    Fly_Struct(i).Ang_veluf=Ang_velUF;
    Fly_Struct(i).Ang_vel=Fil_AV;
end
plot(Fly_Struct(1).Ang_vel)
figure
plot(Fly_Struct(1).Ang_veluf)
% section for box plots included in this function


for j=1:length(Fly_Struct)
    mean_ang_vel(j)=mean(abs(Fly_Struct(j).Ang_vel));
end

figure
boxplot(mean_ang_vel)
xlabel('7.5')
ylabel('Angular Vel')
title('Mean Angular Velocity for Paper Patter 7.5 deg')

end