git remote add private-remote git@github.com:USER_NAME/PRIVATE_REPO.git

function varargout = MotuSense360GUI(varargin)
% MOTUSENSE360GUI MATLAB code for MotuSense360GUI.fig
%      MOTUSENSE360GUI, by itself, creates a new MOTUSENSE360GUI or raises the existing
%      singleton*.
%
%      H = MOTUSENSE360GUI returns the handle to a new MOTUSENSE360GUI or the handle to
%      the existing singleton*.
%
%      MOTUSENSE360GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTUSENSE360GUI.M with the given input arguments.
%
%      MOTUSENSE360GUI('Property','Value',...) creates a new MOTUSENSE360GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MotuSense360GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MotuSense360GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MotuSense360GUI

% Last Modified by GUIDE v2.5 10-Apr-2016 14:57:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MotuSense360GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MotuSense360GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before MotuSense360GUI is made visible.
function MotuSense360GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MotuSense360GUI (see VARARGIN)

% Choose default command line output for MotuSense360GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes MotuSense360GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = MotuSense360GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.txt','Text Files'});

if FileName~=0;

fid = fopen([PathName FileName]);
if fid ==-1
    disp('File open not successful')
    disp(FileName)
    
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
        if isempty(markertimes) ~=1
         elapsed=elapsed-(count*.2); %gets rid of extra time added by the event marking feature (200ms delay)
        end
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


%save data to handles
handles.x=x;
handles.y=y;
handles.z=z;
handles.yaw=yaw;
handles.pitch=pitch;
handles.roll=roll;
handles.t=t;
handles.markertimes=markertimes;
handles.samplingrate=samplingrate;
handles.FileName=FileName;
set(handles.slider1,'Min',t(1));
set(handles.slider1,'Max',t(end));
set(handles.slider1,'Value',t(2));
set(handles.Mintext,'String',['Min:',num2str(t(1))])
set(handles.Maxtext,'String',['Max:',num2str(t(end))])


% samplePeriod=handles.samplingrate;
time=handles.t;
gyrX=handles.yaw;
gyrY = handles.pitch;
gyrZ = handles.roll;
accX = handles.x;
accY = handles.y;
accZ = handles.z;
stopTime=elapsed;


%save changes made to handles
guidata(handles.figure1,handles) 

% samplingrate
% elapsed
% number_of_samples=length(x)
% yes_no=isequal(handles.x,x)
% 
% handles

%play 3D plots
% velocities
end





% --- Executes on button press in AxisPlots.
function AxisPlots_Callback(hObject, eventdata, handles)
% hObject    handle to AxisPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x=handles.x;
y=handles.y;
z=handles.z;
yaw=handles.yaw;
pitch=handles.pitch;
roll=handles.roll;
samplingrate=handles.samplingrate;
markertimes=handles.markertimes;
t=handles.t;

disp([handles.FileName,'Axis Plots'])
% f=figure;
% set(f,'name',[handles.FileName,'Axis Plots'],'numbertitle','off');

% x------------------------------------------------------------------------

%subplot(3,2,1)
axes(handles.axes2);
cla
plot(t,x,'r');
title('X Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('G-value (m^2/s)','Fontsize',15)
axis(handles.axes2,[t(1) t(end) min(x)-1 max(x)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i);
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[xpks,xlocs]=findpeaks(x,'MINPEAKDISTANCE',100);
% hold on
% plot(,t(xlocs),x(xlocs),'ko')

% y------------------------------------------------------------------------

%subplot(3,2,3)
axes(handles.axes3);
cla
plot(t,y,'g');
title('Y Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('G-value (m^2/s)','Fontsize',15)
axis(handles.axes3,[t(1) t(end) min(y)-1 max(y)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[ypks,ylocs]=findpeaks(y,'MINPEAKDISTANCE',100);
% hold on
% plot(t(ylocs),y(ylocs),'ko')

% z------------------------------------------------------------------------

%subplot(3,2,5)
axes(handles.axes4);
cla
plot(t,z,'b');
title('Z Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('G-value (m^2/s)','Fontsize',15)
axis(handles.axes4,[t(1) t(end) min(z)-1 max(z)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[zpks,zlocs]=findpeaks(z,'MINPEAKDISTANCE',100);
% hold on
% plot(t(zlocs),z(zlocs),'ko')


% yaw----------------------------------------------------------------------

%subplot(3,2,2)
axes(handles.axes5);
cla
plot(t,yaw,'r');
title('Yaw Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Rotation (deg/sec)','Fontsize',15)
axis(handles.axes5,[t(1) t(end) min(yaw)-1 max(yaw)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[yawpks,yawlocs]=findpeaks(yaw,'MINPEAKDISTANCE',100);
% hold on
% plot(t(yawlocs),yaw(yawlocs),'ko')

% pitch--------------------------------------------------------------------

%subplot(3,2,4)
axes(handles.axes6);
cla
plot(t,pitch,'g');
title('Pitch Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Rotation (deg/sec)','Fontsize',15)
axis(handles.axes6,[t(1) t(end) min(pitch)-1 max(pitch)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i);  
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[pitchpks,pitchlocs]=findpeaks(pitch,'MINPEAKDISTANCE',100);
% hold on
% plot(t(pitchlocs),pitch(pitchlocs),'ko')

% roll---------------------------------------------------------------------

%subplot(3,2,6)
axes(handles.axes7);
cla
plot(t,roll,'b');
title('Roll Axis','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Rotation (deg/sec)','Fontsize',15)
axis(handles.axes7,[t(1) t(end) min(roll)-1 max(roll)+1])

if isempty(markertimes)==0
for i=1:length(markertimes)
    hold on
    marker=markertimes(i); 
    plot([marker marker],get(gca,'ylim'),'Color','c')
end
end

[rollpks,rolllocs]=findpeaks(roll,'MINPEAKDISTANCE',100);
% hold on
% plot(t(rolllocs),roll(rolllocs),'ko')




% --- Executes on button press in SpectrogramPlots.
function SpectrogramPlots_Callback(hObject, eventdata, handles)
% hObject    handle to SpectrogramPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x=handles.x;
y=handles.y;
z=handles.z;
yaw=handles.yaw;
pitch=handles.pitch;
roll=handles.roll;
t=handles.t;
samplingrate=handles.samplingrate;

% f=figure;
% set(f,'name',[handles.FileName,'Spectrogram Plots'],'numbertitle','off');

% xSpectrograph---------------------------------------------------------------------
% plotting of the spectrogram
%subplot(3,2,1)
axes(handles.axes2);
cla
spectrogram(x, 32, [], [], samplingrate, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude (dB)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency (Hz)')
title('X Axis Spectrogram')


% ySpectrograph---------------------------------------------------------------------


% plotting of the spectrogram
%subplot(3,2,2)
axes(handles.axes3);
cla
spectrogram(y, 32, [], [], samplingrate, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude (dB)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency (Hz)')
title('Y Axis Spectrogram')


% zSpectrograph---------------------------------------------------------------------

% plotting of the spectrogram
%subplot(3,2,3)
axes(handles.axes4);
cla
spectrogram(z, 32, [], [], samplingrate, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude (dB)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency (Hz)')
title('Z Axis Spectrogram')


% yawSpectrograph---------------------------------------------------------------------

% plotting of the spectrogram
%subplot(3,2,4)
axes(handles.axes5);
cla
spectrogram(yaw, 32, [], [], samplingrate, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude (dB)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency (Hz)')
title('Yaw Axis Spectrogram')


% pitchSpectrograph---------------------------------------------------------------------

% plotting of the spectrogram
%subplot(3,2,5)
axes(handles.axes6);
cla
spectrogram(pitch, 32, [], [], samplingrate, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude (dB)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency (Hz)')
title('Pitch Axis Spectrogram')


% rollSpectrograph--------------------------------------------------------

% plotting of the spectrogram
%subplot(3,2,6)
axes(handles.axes7);
cla
spectrogram(roll, 32, [], [], samplingrate, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude (dB)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency (Hz)')
title('Roll Axis Spectrogram')



% --- Executes on button press in StereotopyClassification.
function StereotopyClassification_Callback(hObject, eventdata, handles)
% hObject    handle to StereotopyClassification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x=handles.x;
y=handles.y;
z=handles.z;
yaw=handles.yaw;
pitch=handles.pitch;
roll=handles.roll;
t=handles.t;
samplingrate=handles.samplingrate;

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





xampmeans=[];%preallocate array to store the mean of the peak values from each window
yampmeans=[];
zampmeans=[];
yawampmeans=[];
pitchampmeans=[];
rollampmeans=[];

for i=1:c
xampmeans(i)=mean(xpksstruct{i});
yampmeans(i)=mean(ypksstruct{i});
zampmeans(i)=mean(zpksstruct{i});
yawampmeans(i)=mean(yawpksstruct{i});
pitchampmeans(i)=mean(pitchpksstruct{i});
rollampmeans(i)=mean(rollpksstruct{i});
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




axes(handles.axes8)
horzline=ones(1,length(t))*5;
plot(t,horzline,'k')
title('Stereotopy Classifications','Fontsize',21)
xlabel('Time (s)','Fontsize',18)
set(gca,'YTick',[])
set(gca,'YTickLabel',[])
axis(handles.axes8,[0 t(end)+1 min(horzline)-1 max(horzline)+1])


stereocount=0;

%frame1-frame2 thresholds are used to compare the axis' amplitudes or frequencies to baselines and also to other sterotopy types' avg amp & freq. 
%               The threshold is the std of the amp or freq for a particular axis and particular stereotopy type.

%Wringing------------------------------------------------------------------
% Wringing - small, fast oscillations on x, y, z; high amplitude oscillations on yaw; some oscillations on pitch and roll
% lime green [0.83,0.98,0]

    for i=2:length(yawampmeans)
        if yawampmeans(i)<737 &&yawampmeans(i)> 414.... %|| (yawampmeans(i)-yawampmeans(i-1))>63....
                && yfreqvec(i)<1.5 && yfreqvec(i)>0.4&& xfreqvec(i)<1.6 && xfreqvec(i)>0.5 && zfreqvec(i)< 1.6 && zfreqvec(i)>.5% && fast oscillations xyz
            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
           
           line(tyaw,horzlinenew,'LineWidth',10, 'Color',[0.83,0.98,0]);
           horzlinenew(1)=5.3;
           text(tyaw(1),horzlinenew(1),'Wringing','Fontsize',12);
           
          string=sprintf('Wringing occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
           
stereocount=stereocount+1;
        end
        
        
    end


% %Clapping------------------------------------------------------------------
% % Clapping - fast oscillations on y axis; high amplitude on y, z, and yaw axis?
% %purple blue [0.47,0.46,0.74]

    for i=2:length(yawampmeans)
%         if (yawampmeans(i)-yawampmeans(i-1))>50 && (yampmeans(i)-yampmeans(i-1))>5....
%                 && (zampmeans(i)-zampmeans(i-1))>5;% && fast oscillations on all axis
        if yfreqvec(i)<1.6 && yfreqvec(i)>0.4 && pitchfreqvec(i)<1.8 && pitchfreqvec(i)> 0.60 && yawampmeans(i)< 260 && yawampmeans(i) > 60....
                 && zfreqvec(i)<1.42 && zfreqvec(i)>1.3
%             yampmeans(i)<1.5 &&yampmeans(i)>0.4...
%                 && zampmeans(i)<2.7 && zampmeans(i)>1....
%               && xfreqvec(i)<1.4 && xfreqvec(i)> 0.99&& zfreqvec(i)<1.3 && zfreqvec(i)>1.02...
%             && yawfreqvec(i)<1.3 && yawfreqvec(i)>0.7 && pitchfreqvec(i)<1.4 && pitchfreqvec(i)> 0.87&& rollfreqvec(i)<1.33 && rollfreqvec(i)>0.94
             
            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
            line(tyaw,horzlinenew,'LineWidth',10, 'Color',[0.47,0.46,0.74]);
            horzlinenew(1)=5.3;
            text(tyaw(1),horzlinenew(1),'Clapping','Fontsize',12);
            
          string=sprintf('Clapping occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
          
          
stereocount=stereocount+1;

          
        end
    end
         

% %Tapping-------------------------------------------------------------------
% % Tapping - minimal movement on all axis
% %pink  [0.97,0,0.701]

    for i=2:length(yawampmeans)
%         if (xampmeans(i)-xampmeans(i-1))<10 && (yampmeans(i)-yampmeans(i-1))<10....
%             && (zampmeans(i)-zampmeans(i-1))<10 (yawampmeans(i)-yawampmeans(i-1))<10....
%             &&(pitchampmeans(i)-pitchampmeans(i-1))<10 && (rollampmeans(i)-rollampmeans(i-1))<10;%

if    yawampmeans(i)<72 && yawampmeans(i)>8 && pitchampmeans(i)<100 && pitchampmeans(i)>-35
%     rollampmeans(i)<30 && rollampmeans(i)>5
%     
%     yampmeans(i)<0.5 && yampmeans(i)>0.01 && xampmeans(i)<0.3 && xampmeans(i)>0.05 && zampmeans(i)<0.5 && zampmeans(i)>0.1...
%             && yawampmeans(i)<72 && yawampmeans(i)>8 && pitchampmeans(i)<100 && pitchampmeans(i)>-35 && rollampmeans(i)<30 && rollampmeans(i)>5;        

            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
            line(tyaw,horzlinenew,'LineWidth',10, 'Color',[0.97,0,0.701]);
            horzlinenew(1)=5.3;
            text(tyaw(1),horzlinenew(1),'Tapping','Fontsize',12);
            
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
%         if (zampmeans(i)-zampmeans(i-1))>50 && (pitchampmeans(i)-pitchampmeans(i-1))>50; % && fast oscillations on all axis

    if   zampmeans(i)<20 && zampmeans(i)>1.8 &&pitchampmeans(i)<550 &&pitchampmeans(i)>140
%     && yfreqvec(i)< && yfreqvec(i)> && xfreqvec(i)< && xfreqvec(i)> && zfreqvec(i)< && zfreqvec(i)>...
%             && yawfreqvec(i)< && yawfreqvec(i)> && pitchfreqvec(i)< && pitchfreqvec(i)> && rollfreqvec(i)< && rollfreqvec(i)>

            tyaw=t(yawpkslocsstruct{i});
            horzlinenew=horzline(yawpkslocsstruct{i});
            line(tyaw,horzlinenew,'LineWidth',10, 'Color',[.6,0.3,0.3]);
            horzlinenew(1)=5.3;
            if stereocount>1
                horzlinenew(1)=4.6;
            end
            text(tyaw(1),horzlinenew(1),'Flapping','Fontsize',12);
            
          string=sprintf('Flapping occured from %.2f seconds to %.2f seconds',tyaw(1),tyaw(end));
          disp(string)
          
          stereocount=stereocount+1;


        end
    end


function Window_Callback(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Window as text
%        str2double(get(hObject,'String')) returns contents of Window as a double
handles.windowval=str2double(get(hObject,'String'));
guidata(handles.figure1,handles) 



% --- Executes during object creation, after setting all properties.
function Window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Apply.
function Apply_Callback(hObject, eventdata, handles)
% hObject    handle to Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.windowval;
minval=get(handles.slider1,'Value');
maxval=get(handles.slider1,'Value')+handles.windowval;

[~, startloc] = min(abs(handles.t - minval));
[~, endloc] = min(abs(handles.t - maxval));

handles.t=handles.t(startloc:endloc);
handles.x=handles.x(startloc:endloc);
handles.y=handles.y(startloc:endloc);
handles.z=handles.z(startloc:endloc);
handles.yaw=handles.yaw(startloc:endloc);
handles.pitch=handles.pitch(startloc:endloc);
handles.roll=handles.roll(startloc:endloc);

AxisPlots_Callback(hObject,eventdata,handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderval=get(hObject,'Value');
set(handles.Valuetext,'String',['Value:',num2str(sliderval)]);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
