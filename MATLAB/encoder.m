function tick_to_meter = encoder(tick)
%This function takes in the encoder tick reading(s) then converts them into
%length [m]. The input column vecors: tick(:,1) is time, tick(:,2) is
%the x1 cart tick position, and tick(:,3) is the x2 cart tick position. The
%output tick_to_meter gives back the time vector plus the converted tick
%vectors.

conv = 1/4096; %[rev/pulses]
circ = 0.1; %circumference [m]

tick_to_meter(:,1) = tick(:,1); %time [s]
tick_to_meter(:,2:3) = circ*conv*tick(:,2:3); %ticks (or pulses) to m [m/pulse]
end