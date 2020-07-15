%% Sleep Teensy sandbox:

addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2\Control_scripts\Sleep_state'))

%% initialize
ports= seriallist;

teensy = serial(ports{end}, 'Baudrate', 9600);
fopen(teensy)
flushinput(teensy)
fprintf('Now connected to %s \n', ports{end});


%% another attempt

coms = instrfind;

if ~isempty(coms)
    fclose(coms);
end

SerialID = serial('COM3', 'BaudRate', 19200, 'InputBufferSize', 12000);
fopen(SerialID);

figure(1)
Fs = 2000; 
eT = 0;
total_T = 100; 
plot((eT:(1/Fs):total_T),nan(1,length(eT:(1/Fs):total_T)))

while eT <= total_T
    this_T = eT; 
    eT = eT+1; 
    tvec = this_T:(1/Fs):eT;
    tvec = tvec(1:end-1);
this_data = fread(SerialID,2000);
plot(tvec, this_data');
drawnow
pause(2)
end

%%
serialportlist("available");
arduinoObj = serialport("COM3",9600);
configureTerminator(arduinoObj,"CR/LF");
flush(arduinoObj);
arduinoObj.UserData = struct("Data",[],"Count",1);
%%
figure(2)
close(2)
figure(2)
Fs = 2000; 
eT = 1;
print_t = 10; % print every this many samples
total_T = 10000; 
this_data = []; tvec = []; 
% plot((eT:(1/Fs):total_T),nan(1,length(eT:(1/Fs):total_T)))
count = 0; 
flush(arduinoObj);

h = animatedline('MaximumNumPoints',1000);
axis([0,1,-1,1])

while eT < total_T
    count = count+1; 
    tvec(count,1) = eT/1000;
    this_data(count,1) = str2double(readline(arduinoObj));
%     if abs(tvec(1)- tvec(end)) == print_t
        fprintf('%0.1f - %0.1f - %0.1f\n', eT, tvec(end), this_data(end))
        addpoints(h, tvec, this_data)
%             plot(tvec, this_data)
        drawnow
        tvec = [];
        this_data = [];
%     pause(0.005)
% %     hold on
%     drawnow
%     end
        eT = eT+1;

end

%%
while eT <= total_T
    this_T = eT; 
    eT = eT+1; 
    tvec = this_T:(1/Fs):eT;
    tvec = tvec(1:end-1);
this_data = str2double(readline(arduinoObj));
plot(tvec, this_data');
drawnow
pause(2)
end


% readTeensyData(arduinoObj)
for ii = 1:3
    data
% configureCallback(arduinoObj,"terminator",@readTeensyData);
hold on
pause(5)
ii = ii +1;
end

%%

s = serialport("COM3",9600);

for ii = 1:10
configureCallback(s,"byte",100,@readTeensyData)


end

configureCallback(s,"off")

