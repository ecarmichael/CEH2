%%% control script for tracking wheel position and giving rewards using
%%% neuralynx






%% init

max_time = 10; % time in seconds; 
max_laps = 50; % max number of laps

% port names
pulse_cmd = '-DigitalIOTtlPulse AcqSystem1_0 '; 
TTL_IO_port = 0; 
TTL_lick = 0;
TTL_clk_a = 2;
TTL_clk_b = 3;



%tone
amp = 10; fs = 20500; duration = 0.25; freq = 2000; values = 0:1/fs:duration; 
a = amp*sin(2*pi*freq*values);



%% Connect to Cheetah
connected = NlxAreWeConnected();
if connected ==1
    disp 'We are already connected';
else
    serverName = '192.168.3.100';
    disp(sprintf('Connecting to %s...', serverName));
    success = NlxConnectToServer(serverName);
    if success ~= 1
        disp(sprintf('FAILED connect to %s.', serverName));
        return;
    else
        disp(sprintf('Connected to Cheetah at %s.', serverName));
    end
end
%% Populate fields

TotalLaps = 0;
CorrectLaps = 0;
IncorrectLaps = 0;
state = 0;

lick = 0;


%% make the tracking plot
f1_handle = figure(1);
clf
fs = 18;
set(f1_handle,'Position',[400,100,900,900])
t_handle = title('Task starting...'); set(t_handle,'FontSize',fs);


% trial counts
laps_handle = text(0.2,0.9,sprintf('Correct Laps %d | Total Laps: %d',TotalLaps,CorrectLaps),'Interpreter','none');
set(laps_handle,'FontSize',fs);

subplot(2,2,1)
cla
axis off

% make the circle that will rotate
subplot(2,2,2)
cla
axis off
hold on
rho = ones(24,1); 
step = 24; 
theta = 0:360/step:(360 - 360/step);
x = [0 rho(1)]; 
y = [0 0];

for ii = 1:step
    y2=0+(step*sind(theta(ii)));
    x2=0+(step*cosd(theta(ii)));
    plot([0 x2],[0 y2], 'k')
    text(
    
end


v = version('-release'); 




if str2double(v(1:4)) < 2025
    [x,y] = pol2cart(theta * pi/180, rho);
    c = compass(x, y, 'k');
    for ii = 1:step
        if ii ==1 
        c(ii).LineWidth = 2; 
        c(ii).Color = "Red"; 
        c(ii).Marker = 'o';
        else
        c(ii).Marker = '.'; 
        c(ii).Color = "Black"; 
        end
    end
else
    compassplot(theta, rho);
end

subplot(2,2,3:4)
cla
ylabel('lick rate')
xlabel('time')


%% run the program
t0 = clock; 
first = 1; 
while (etime(clock, t0) < max_time) && (TotalLaps  < max_laps)

    % disp(etime(clock, t0))
    
    [evt] = NLX_IO_in(TTL_IO_port);
    
    if evt == 1
        lick = 1; 
    end

    if lick ==1
        if l_on == 0 
            fprintf('Lick at t: %0.2f', etime(clock, t0))
            l_on = 1; first = 0; 
        end
    elseif ~ evt
        if l_on == 1 %&& first ~= 1
        fprintf('Lick off at t: %0.2f\n', etime(clock, t0))
        l_on = 0;
        end
    end
        
end
