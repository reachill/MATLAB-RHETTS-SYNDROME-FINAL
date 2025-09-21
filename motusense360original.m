%preloaded with Nicoles data


close all; clear all; clc
fid = fopen('/Users/ntunick19/Desktop/nomo_wringing_tapping_clapping_4_14.txt'); 
if fid ==-1
    disp('File open not successful')
    
else
    disp('File open successful')
    while feof(fid) == 0

        % Read one line into a string variable
        [A] = fscanf(fid,'%f'); %every six numbers are x,y,z,yaw,pitch,roll
        count=A(end)-1;
        A(end)=[];
        markertimes=A(end-count:end)/1000;
        A(end-count:end)=[];
        elapsed=A(end)/1000;
        elapsed=elapsed-(count*.2); %gets rid of extra time added by the event marking feature (200ms delay)
        A(end)=[];
        
        while rem(length(A),6) ~=0
            A(1)=[];
        end

        acclgyro = reshape(A,6,length(A)/6);
    end
        closeresult= fclose(fid);

    
    if closeresult ==0
        disp('File close successful')
    else
        disp('File close not succesful')
    end
end

% zerorange=t(0:


samplingrate=length(acclgyro(1,:))/elapsed;
t=1/(samplingrate):1/(samplingrate):length(acclgyro(1,:))/(samplingrate);

%accelerometer changes-----------------------------------------------------
x= acclgyro(1,:)/256;
y=acclgyro(2,:)/256;
z=(acclgyro(3,:)/256)-1;
accl=[x,y,z];

%gyroscope changes---------------------------------------------------------
yaw=acclgyro(4,:)/14.375;
pitch= acclgyro(5,:)/14.375;
roll= acclgyro (6,:)/14.375;
gyro=[yaw,pitch,roll];

%%
% x------------------------------------------------------------------------

subplot(3,2,1)
plot(t,x,'r');
title('X Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('G-value (m^2/s)','Fontsize',15)
axis([0 t(end) min(x)-1 max(x)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i);
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[xpks,xlocs]=findpeaks(x,'MINPEAKDISTANCE',100);
hold on
% plot(t(xlocs),x(xlocs),'ko')

% y------------------------------------------------------------------------

subplot(3,2,3)
plot(t,y,'g');
title('Y Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('G-value (m^2/s)','Fontsize',15)
axis([0 t(end) min(y)-1 max(y)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[ypks,ylocs]=findpeaks(y,'MINPEAKDISTANCE',100);
hold on
% plot(t(ylocs),y(ylocs),'ko')

% z------------------------------------------------------------------------

subplot(3,2,5)
plot(t,z,'b');
title('Z Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('G-value (m^2/s)','Fontsize',15)
axis([0 t(end) min(z)-1 max(z)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[zpks,zlocs]=findpeaks(z,'MINPEAKDISTANCE',100);
hold on
% plot(t(zlocs),z(zlocs),'ko')


% yaw----------------------------------------------------------------------

subplot(3,2,2)
plot(t,yaw,'r');
title('Yaw Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Rotation (deg/sec)','Fontsize',15)
axis([0 t(end) min(yaw)-1 max(yaw)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[yawpks,yawlocs]=findpeaks(yaw,'MINPEAKDISTANCE',100);
hold on
% plot(t(yawlocs),yaw(yawlocs),'ko')

% pitch--------------------------------------------------------------------

subplot(3,2,4)
plot(t,pitch,'g');
title('Pitch Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Rotation (deg/sec)','Fontsize',15)
axis([0 t(end) min(pitch)-1 max(pitch)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i);  
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[pitchpks,pitchlocs]=findpeaks(pitch,'MINPEAKDISTANCE',100);
hold on
% plot(t(pitchlocs),pitch(pitchlocs),'ko')

% roll---------------------------------------------------------------------

subplot(3,2,6)
plot(t,roll,'b');
title('Roll Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Rotation (deg/sec)','Fontsize',15)
axis([0 t(end) min(roll)-1 max(roll)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[rollpks,rolllocs]=findpeaks(roll,'MINPEAKDISTANCE',100);
hold on
% plot(t(rolllocs),roll(rolllocs),'ko')
set(gcf,'color','w');
% 
% % xSpectrograph---------------------------------------------------------------------
% 
% 
% 
% % time & discretisation parameters
% N = length(x);       
% 
% % plotting of the waveform
% figure
% subplot(2,1,1)
% plot(t, x, 'r')
% xlim([0 max(t)])
% ylim([-1.1*max(abs(x)) 1.1*max(abs(x))])
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Normalized amplitude')
% title('The signal in the time domain')
% 
% % plotting of the spectrogram
% subplot(2,1,2)
% spectrogram(x, 32, [], [], samplingrate, 'yaxis')
% h = colorbar;
% set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(h, 'Magnitude, dB');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Spectrogram of the signal')
% set(gcf,'color','w');
% % ySpectrograph---------------------------------------------------------------------
% 
% 
% 
% % time & discretisation parameters
% N = length(y);       
% 
% % plotting of the waveform
% figure
% subplot(2,1,1)
% plot(t, y, 'r')
% xlim([0 max(t)])
% ylim([-1.1*max(abs(y)) 1.1*max(abs(y))])
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Normalized amplitude')
% title('The signal in the time domain')
% 
% % plotting of the spectrogram
% subplot(2,1,2)
% spectrogram(y, 32, [], [], samplingrate, 'yaxis')
% h = colorbar;
% set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(h, 'Magnitude, dB');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Spectrogram of the signal')
% 
% 
% set(gcf,'color','w');
% % zSpectrograph---------------------------------------------------------------------
% 
% 
% 
% % time & discretisation parameters
% N = length(z);       
% 
% % plotting of the waveform
% figure
% subplot(2,1,1)
% plot(t, z, 'r')
% xlim([0 max(t)])
% ylim([-1.1*max(abs(z)) 1.1*max(abs(z))])
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Normalized amplitude')
% title('The signal in the time domain')
% 
% % plotting of the spectrogram
% subplot(2,1,2)
% spectrogram(z, 32, [], [], samplingrate, 'yaxis')
% h = colorbar;
% set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(h, 'Magnitude, dB');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Spectrogram of the signal')
% 
% 
% set(gcf,'color','w');
% % yawSpectrograph---------------------------------------------------------------------
% 
% 
% 
% % time & discretisation parameters
% N = length(yaw);       
% 
% % plotting of the waveform
% figure
% subplot(2,1,1)
% plot(t, yaw, 'r')
% xlim([0 max(t)])
% ylim([-1.1*max(abs(yaw)) 1.1*max(abs(yaw))])
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Normalized amplitude')
% title('The signal in the time domain')
% 
% % plotting of the spectrogram
% subplot(2,1,2)
% spectrogram(yaw, 32, [], [], samplingrate, 'yaxis')
% h = colorbar;
% set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(h, 'Magnitude, dB');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Spectrogram of the signal')
% set(gcf,'color','w');
% 
% % pitchSpectrograph---------------------------------------------------------------------
% 
% 
% 
% % time & discretisation parameters
% N = length(pitch);       
% 
% % plotting of the waveform
% figure
% subplot(2,1,1)
% plot(t, pitch, 'r')
% xlim([0 max(t)])
% ylim([-1.1*max(abs(pitch)) 1.1*max(abs(pitch))])
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Normalized amplitude')
% title('The signal in the time domain')
% 
% % plotting of the spectrogram
% subplot(2,1,2)
% spectrogram(pitch, 32, [], [], samplingrate, 'yaxis')
% h = colorbar;
% set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(h, 'Magnitude, dB');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Spectrogram of the signal')
% set(gcf,'color','w');
% 
% % rollSpectrograph---------------------------------------------------------
% 
% 
% 
% % time & discretisation parameters
% N = length(roll);       
% 
% % plotting of the waveform
% figure
% subplot(2,1,1)
% plot(t, roll, 'r')
% xlim([0 max(t)])
% ylim([-1.1*max(abs(roll)) 1.1*max(abs(roll))])
% grid on
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Normalized amplitude')
% title('The signal in the time domain')
% 
% % plotting of the spectrogram
% subplot(2,1,2)
% spectrogram(roll, 32, [], [], samplingrate, 'yaxis')
% h = colorbar;
% set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(h, 'Magnitude, dB');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Spectrogram of the signal')
% set(gcf,'color','w');
%%
%--------------------------------------------------------------------------


%create windows of time
%find std of amplitudes and campare them to previous windows
%find osillations/frequency/ period of amplitudes and compare


% Wringing - small, fast oscillations on x, y, z; high amplitude oscillations on yaw; some oscillations on pitch and roll
% 
% Clapping - fast oscillations on all axis; high amplitude on y, z, and yaw axis?
% 
% Tapping - minimal movement on all axis
% 
% Flapping - fast oscillations on all axis; high amplitude on z and pitch
% 
% Hand-mouthing - change in x and z baseline amplitude; small oscillations in x, medium oscillations in y and z; fast oscillations
%                 in yaw, pitch, and roll, medium sized amplitude in yaw, varying peaks in pitch and roll




%Separate into 5 second windows and find pks in those windows--------------
 %(split columns from frames speraely then find pks)
 
framewidth=820;% splits frames into about 5 seconds worth of data, 820 elements is about 5 seconds
overlap=328;%overlaps about 2 seconds
xframes=buffer(x,framewidth); 
yframes=buffer(y,framewidth);
zframes=buffer(z,framewidth);
yawframes=buffer(yaw,framewidth);
pitchframes=buffer(pitch,framewidth);
rollframes=buffer(roll,framewidth);

[r c]=size(xframes);
xstruct=cell(1,c);
ystruct=cell(1,c);
zstruct=cell(1,c);
yawstruct=cell(1,c);
pitchstruct=cell(1,c);
rollstruct=cell(1,c);

for i=1:c
    xstruct{i}=xframes(:,i);
    ystruct{i}=yframes(:,i);
    zstruct{i}=zframes(:,i);
    yawstruct{i}=yawframes(:,i);
    pitchstruct{i}=pitchframes(:,i);
    rollstruct{i}=rollframes(:,i);
    
end


%amplitudes----------------------------------------------------------------

xpksstruct=cell(1,c); %preallocate structures to store peak values of each window
xpkslocsstruct=cell(1,c);

ypksstruct=cell(1,c);
ypkslocsstruct=cell(1,c);

zpksstruct=cell(1,c);
zpkslocsstruct=cell(1,c);

yawpksstruct=cell(1,c);
yawpkslocsstruct=cell(1,c);

pitchpksstruct=cell(1,c);
pitchpkslocsstruct=cell(1,c);
   
rollpksstruct=cell(1,c);
rollpkslocsstruct=cell(1,c);

for i=1:c
    [xpksstruct{i},xpkslocsstruct{i}]=findpeaks(xstruct{i},'MINPEAKDISTANCE',100);
    [ypksstruct{i},ypkslocsstruct{i}]=findpeaks(ystruct{i},'MINPEAKDISTANCE',100);
    [zpksstruct{i},zpkslocsstruct{i}]=findpeaks(zstruct{i},'MINPEAKDISTANCE',100);
    [yawpksstruct{i},yawpkslocsstruct{i}]=findpeaks(yawstruct{i},'MINPEAKDISTANCE',100);
    [pitchpksstruct{i},pitchpkslocsstruct{i}]=findpeaks(pitchstruct{i},'MINPEAKDISTANCE',100);
    [rollpksstruct{i},rollpkslocsstruct{i}]=findpeaks(rollstruct{i},'MINPEAKDISTANCE',100);  
    
    xpkslocsstruct{i}=xpkslocsstruct{i}+((i-1)*820);%this makes the indexes for each frame sequential to be in the same order as the indexes in t
    ypkslocsstruct{i}=ypkslocsstruct{i}+((i-1)*820);
    zpkslocsstruct{i}=zpkslocsstruct{i}+((i-1)*820);
    yawpkslocsstruct{i}=yawpkslocsstruct{i}+((i-1)*820);
    pitchpkslocsstruct{i}=pitchpkslocsstruct{i}+((i-1)*820);
    rollpkslocsstruct{i}=rollpkslocsstruct{i}+((i-1)*820);
 
    n=(820*c)-length(t);%this section gets rid of the extra zeros added in the frames by the buffer function
    m=(820*c)-n;
    xpkslocsstruct{i}(m:end)=[];
    ypkslocsstruct{i}(m:end)=[];
    zpkslocsstruct{i}(m:end)=[];
    yawpkslocsstruct{i}(m:end)=[];
    pitchpkslocsstruct{i}(m:end)=[];
    rollpkslocsstruct{i}(m:end)=[];
    
end


%--------------------------------------




xpitstruct=cell(1,c); %preallocate structures to store peak values of each window
xpitlocsstruct=cell(1,c);

ypitstruct=cell(1,c);
ypitlocsstruct=cell(1,c);

zpitstruct=cell(1,c);
zpitlocsstruct=cell(1,c);

yawpitstruct=cell(1,c);
yawpitlocsstruct=cell(1,c);

pitchpitstruct=cell(1,c);
pitchpitlocsstruct=cell(1,c);
   
rollpitstruct=cell(1,c);
rollpitlocsstruct=cell(1,c);

for i=1:c
    [xpitstruct{i},xpitlocsstruct{i}]=findpeaks(-1*xstruct{i},'MINPEAKDISTANCE',100);
    [ypitstruct{i},ypitlocsstruct{i}]=findpeaks(-1*ystruct{i},'MINPEAKDISTANCE',100);
    [zpitstruct{i},zpitlocsstruct{i}]=findpeaks(-1*zstruct{i},'MINPEAKDISTANCE',100);
    [yawpitstruct{i},yawpitlocsstruct{i}]=findpeaks(-1*yawstruct{i},'MINPEAKDISTANCE',100);
    [pitchpitstruct{i},pitchpitlocsstruct{i}]=findpeaks(-1*pitchstruct{i},'MINPEAKDISTANCE',100);
    [rollpitstruct{i},rollpitlocsstruct{i}]=findpeaks(-1*rollstruct{i},'MINPEAKDISTANCE',100);  
    
        xpitstruct{i}=-1*xpitstruct{i};
    ypitstruct{i}=-1*ypitstruct{i};
    zpitstruct{i}=-1*zpitstruct{i};
    yawpitstruct{i}=-1*yawpitstruct{i};
    pitchpitstruct{i}=-1*pitchpitstruct{i};
    rollpitstruct{i}=-1*rollpitstruct{i};
    
    xpitlocsstruct{i}=xpitlocsstruct{i}+((i-1)*820);%this makes the indexes for each frame sequential to be in the same order as the indexes in t
    ypitlocsstruct{i}=ypitlocsstruct{i}+((i-1)*820);
    zpitlocsstruct{i}=zpitlocsstruct{i}+((i-1)*820);
    yawpitlocsstruct{i}=yawpitlocsstruct{i}+((i-1)*820);
    pitchpitlocsstruct{i}=pitchpitlocsstruct{i}+((i-1)*820);
    rollpitlocsstruct{i}=rollpitlocsstruct{i}+((i-1)*820);
 
    n=(820*c)-length(t);%this section gets rid of the extra zeros added in the frames by the buffer function
    m=(820*c)-n;
    xpitlocsstruct{i}(m:end)=[];
    ypitlocsstruct{i}(m:end)=[];
    zpitlocsstruct{i}(m:end)=[];
    yawpitlocsstruct{i}(m:end)=[];
    pitchpitlocsstruct{i}(m:end)=[];
    rollpitlocsstruct{i}(m:end)=[];
    
end






%--------------------------------------

xampmeans=[];%preallocate array to store the mean of the peak values from each window
yampmeans=[];
zampmeans=[];
yawampmeans=[];
pitchampmeans=[];
rollampmeans=[];

for i=1:c
    
    while(length(xpksstruct{i})~=length(xpitstruct{i}))
    
    if length(xpksstruct{i}) > length(xpitstruct{i})
        xpksstruct{i}(end)=[];
    end
    
    if length(xpksstruct{i}) < length(xpitstruct{i})
        xpitstruct{i}(end)=[];
    end
    
    end
  
      

        while(length(ypksstruct{i})~=length(ypitstruct{i}))

        if length(ypksstruct{i}) > length(ypitstruct{i})
        ypksstruct{i}(end)=[];
    end
    
    if length(ypksstruct{i}) < length(ypitstruct{i})
        ypitstruct{i}(end)=[];
    end
    
    end
    
    
    
while(length(zpksstruct{i})~=length(zpitstruct{i}))
        if length(zpksstruct{i}) > length(zpitstruct{i})
        zpksstruct{i}(end)=[];
    end
    
    if length(zpksstruct{i}) < length(zpitstruct{i})
        zpitstruct{i}(end)=[];
    end
    
end



while(length(yawpksstruct{i})~=length(yawpitstruct{i}))
if length(yawpksstruct{i}) > length(yawpitstruct{i})
        yawpksstruct{i}(end)=[];
    end
    
    if length(yawpksstruct{i}) < length(yawpitstruct{i})
        yawpitstruct{i}(end)=[];
    end
    end
    
    while(length(pitchpksstruct{i})~=length(pitchpitstruct{i}))
        if length(pitchpksstruct{i}) > length(pitchpitstruct{i})
        pitchpksstruct{i}(end)=[];
    end
    
    if length(pitchpksstruct{i}) < length(pitchpitstruct{i})
        pitchpitstruct{i}(end)=[];
    end
    end
    
   while(length(rollpksstruct{i})~=length(rollpitstruct{i}))
            if length(rollpksstruct{i}) > length(rollpitstruct{i})
        rollpksstruct{i}(end)=[];
    end
    
    if length(rollpksstruct{i}) < length(rollpitstruct{i})
        rollpitstruct{i}(end)=[];
    end
    end
xampmeans(i)=mean(abs(xpksstruct{i}-xpitstruct{i}));
yampmeans(i)=mean(abs(ypksstruct{i}-ypitstruct{i}));
zampmeans(i)=mean(abs(zpksstruct{i}-zpitstruct{i}));
yawampmeans(i)=mean(abs(yawpksstruct{i}-yawpitstruct{i}));
pitchampmeans(i)=mean(abs(pitchpksstruct{i}-pitchpitstruct{i}));
rollampmeans(i)=mean(abs(rollpksstruct{i}-rollpitstruct{i}));
end


%Frequency-----------------------------------------------------------------

%Frequency=1/T
xfreqvec=[];
for i=1:length(xpkslocsstruct)
    xlocsnew=xpkslocsstruct{i};
    xfreqvec(i)=calcfreq(xlocsnew,t);
end

yfreqvec=[];
for i=1:length(ypkslocsstruct)
    ylocsnew=ypkslocsstruct{i};
    
    yfreqvec(i)=calcfreq(ylocsnew,t);
end

zfreqvec=[];
for i=1:length(zpkslocsstruct)
    zlocsnew=zpkslocsstruct{i};
    
    zfreqvec(i)=calcfreq(zlocsnew,t);
end

yawfreqvec=[];
for i=1:length(yawpkslocsstruct)
    yawlocsnew=yawpkslocsstruct{i};
    
    yawfreqvec(i)=calcfreq(yawlocsnew,t);
end

pitchfreqvec=[];
for i=1:length(pitchpkslocsstruct)
    pitchlocsnew=pitchpkslocsstruct{i};
    
    pitchfreqvec(i)=calcfreq(pitchlocsnew,t);
end

rollfreqvec=[];
for i=1:length(rollpkslocsstruct)
    rolllocsnew=rollpkslocsstruct{i};
    
    rollfreqvec(i)=calcfreq(rolllocsnew,t);
end

%%

figure
horzline=ones(1,length(t))*5;
plot(t,horzline,'k')
title('Stereotopy Classifications','Fontsize',22)
xlabel('Time (s)','Fontsize',18)
set(gca,'YTick',[])
set(gca,'YTickLabel',[])
axis([0 t(end)+1 min(horzline)-1 max(horzline)+1])
set(gcf,'color','w');

stereocount=0;

%frame1-frame2 thresholds are used to compare the axis' amplitudes or frequencies to baselines and also to other sterotopy types' avg amp & freq. 
%               The threshold is the std of the amp or freq for a particular axis and particular stereotopy type.

%Wringing------------------------------------------------------------------
% Wringing - small, fast oscillations on x, y, z; high amplitude oscillations on yaw; some oscillations on pitch and roll
% lime green [0.83,0.98,0]

    for i=2:length(yawampmeans)
        if yawampmeans(i)<500 &&yawampmeans(i)> 260.... %|| (yawpksmeans(i)-yawpksmeans(i-1))>63....
                && yfreqvec(i)<1.30 && yfreqvec(i)>0.70&& xfreqvec(i)<1.5 && xfreqvec(i)>0.81 && zfreqvec(i)< 1.6 && zfreqvec(i)>.02% && fast oscillations xyz
            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
           
           line(tyaw,horzlinenew,'LineWidth',50, 'Color',[0.83,0.98,0]);
           horzlinenew(1)=5.3;
           text(tyaw(1),horzlinenew(1),'Wringing','Fontsize',15);
           
          string=sprintf('Wringing occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
           
stereocount=stereocount+1;
        end
        
        
    end


% %Clapping------------------------------------------------------------------
% % Clapping - fast oscillations on y axis; high amplitude on y, z, and yaw axis?
% %purple blue [0.47,0.46,0.74]

    for i=2:length(yawampmeans)
%         if (yawpksmeans(i)-yawpksmeans(i-1))>50 && (ypksmeans(i)-ypksmeans(i-1))>5....
%                 && (zpksmeans(i)-zpksmeans(i-1))>5;% && fast oscillations on all axis
        if yawampmeans(i)< 280 && yawampmeans(i)>60 && yampmeans(i)<2.83 &&yampmeans(i)>-0.5...
                && zampmeans(i)<4.3 && zampmeans(i)>-0.5....
                  && yfreqvec(i)<1.2 && yfreqvec(i)>0.82 && pitchfreqvec(i)<1.4 && pitchfreqvec(i)> 0.87
%               && xfreqvec(i)<1.4 && xfreqvec(i)> 0.99&& zfreqvec(i)<1.3 && zfreqvec(i)>1.02...
%             && yawfreqvec(i)<1.3 && yawfreqvec(i)>0.7 && pitchfreqvec(i)<1.4 && pitchfreqvec(i)> 0.87&& rollfreqvec(i)<1.33 && rollfreqvec(i)>0.94
             
            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
            line(tyaw,horzlinenew,'LineWidth',50, 'Color',[0.47,0.46,0.74]);
            horzlinenew(1)=5.3;
            text(tyaw(1),horzlinenew(1),'Clapping','Fontsize',15);
            
          string=sprintf('Clapping occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
          
          
stereocount=stereocount+1;

          
        end
    end
         

% %Tapping-------------------------------------------------------------------
% % Tapping - minimal movement on all axis
% %pink  [0.97,0,0.701]

    for i=2:length(yawampmeans)
%         if (xpksmeans(i)-xpksmeans(i-1))<10 && (ypksmeans(i)-ypksmeans(i-1))<10....
%             && (zpksmeans(i)-zpksmeans(i-1))<10 (yawpksmeans(i)-yawpksmeans(i-1))<10....
%             &&(pitchpksmeans(i)-pitchpksmeans(i-1))<10 && (rollpksmeans(i)-rollpksmeans(i-1))<10;%

if    yampmeans(i)<1.00 && yampmeans(i)>-0.15 && xampmeans(i)<0.75 && xampmeans(i)>-0.11 && zampmeans(i)<0.14 && zampmeans(i)>-0.3...
            && yawampmeans(i)<66.3 && yawampmeans(i)>13.9 && pitchampmeans(i)<100 && pitchampmeans(i)>10 && rollampmeans(i)<130 && rollampmeans(i)>10;        

            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
            line(tyaw,horzlinenew,'LineWidth',50, 'Color',[0.97,0,0.701]);
            horzlinenew(1)=5.3;
            text(tyaw(1),horzlinenew(1),'Tapping','Fontsize',15);
            
          string=sprintf('Tapping occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
          
          stereocount=stereocount+1;


        end
    end
    
    
% 
% %Flapping------------------------------------------------------------------
% % Flapping - fast oscillations on all axis; high amplitude on z and pitch
% %brown [.6,0.3,0.3]
% 

    for i=2:length(yawampmeans)
%         if (zpksmeans(i)-zpksmeans(i-1))>50 && (pitchpksmeans(i)-pitchpksmeans(i-1))>50; % && fast oscillations on all axis

    if   zampmeans(i)<3 && zampmeans(i)>-0.5 &&pitchampmeans(i)<230 &&pitchampmeans(i)>90
%     && yfreqvec(i)< && yfreqvec(i)> && xfreqvec(i)< && xfreqvec(i)> && zfreqvec(i)< && zfreqvec(i)>...
%             && yawfreqvec(i)< && yawfreqvec(i)> && pitchfreqvec(i)< && pitchfreqvec(i)> && rollfreqvec(i)< && rollfreqvec(i)>

            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
            line(tyaw,horzlinenew,'LineWidth',50, 'Color',[.6,0.3,0.3]);
            horzlinenew(1)=5.3;
%             if stereocount>1
%                 horzlinenew(1)=4.6;
%             end
            text(tyaw(1),horzlinenew(1),'Flapping','Fontsize',15);
            
          string=sprintf('Flapping occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
          
          stereocount=stereocount+1;


        end
    end


% 
% %Hand-Mouthing-------------------------------------------------------------
% % Hand-mouthing - change in x and z baseline amplitude; small oscillations in x, medium oscillations in y and z; fast oscillations
% %                 in yaw, pitch, and roll, medium sized amplitude in yaw, varying peaks in pitch and roll
% %gery [0.5,0.5,0.5]

%     for i=2:length(yawpksmeans)
%         if (xpksmeans(i)-xpksmeans(i-1))>50 && (zpksmeans(i)-zpksmeans(i-1))>50....
%                 && (yawpksmeans(i)-yawpksmeans(i-1))>25;% small osxillations x, med y and z, large yaw pitch roll
%             tyaw=t(yawlocsstruct{i});
%             horzlinenew=horzline(yawlocsstruct{i});
%             line(tyaw,horzlinenew,'LineWidth',3, 'Color',[0.5,0.5,0.5]);
%             horzlinenew(1)=4.9;
%             text(tyaw(1),horzlinenew(1),'Hand\_Mouthing','Fontsize',12);
%             
%           string=sprintf('Hand Mouthing occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
%           disp(string)
%             
%         end
%     end
% 
% 
