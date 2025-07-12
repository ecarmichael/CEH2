
a = arduino("COM4", "UnoR4Minima" );

encoder = rotaryEncoder(a, "D2", "D3", 24)
%% 
d = [];
 brightness_step = (5-0)/100; 
   for ii = 1:100
      writePWMVoltage(a, 'D11', ii*brightness_step);
      d(ii) = readVoltage(a, 'A0');
      pause(0.01);
   end
 
   for ii = 1:20
      writePWMVoltage(a, 'D11', 5-ii*brightness_step);
      pause(0.1);
   end

   figure(10)
   clf
   hold on
   plot(d)
    
   %% read the encoder?
   % configurePin(a, 'D2', 'DigitalInput'); 
   % configurePin(a, 'D3', 'DigitalInput'); 
figure(10)
clf
hold on
c = NaN(1,100);
t = c;
s = c;
e_time = 0; 
tic
while e_time < 10
    e_time = toc;
    plot(e_time, readDigitalPin(a, 'D3'), 'd');

    % plot(e_time, readCount(encoder), 'd');
    % s(ii) = readSpeed(encoder);
end