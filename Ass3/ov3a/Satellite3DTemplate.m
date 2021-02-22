clear all
close all
clc
% 
% pos  = sym('pa',[3,1],'real');
% R = sym('ra',[3,3],'real');
% vel  = sym('va',[3,1],'real');
% ang_vel = sym('aa',[3,1],'real');
% 
% r_o = 6356*10^3 + 36*10^6;
% 
% pos = [r_o; 0; 0];
% 
% 
% 
% 
% % Define your initial state, e.g. as:
% state = [pos;
%          reshape(R,9,1);
%          vel;
%          ang_vel];
%      
% %      parameters = [M;
% %                  m_T;]


%%

r = 6356*10^3 + 36*10^6;
length_sat = 0.5;
density = 1;
m_sat = (0.5)^3*density;

% Define your initial state, e.g. as:
% state = [position;
%          orientation;
%          velocity;
%          angular velocity];

position = [0; 0; r];
orientation = eye(3);
velocity = [0; 0; 1];
omega = [0; 0; 1];

state = [position; reshape(orientation, 9,1); velocity; omega];

%%





l = 36 * 10^6;
M = 6 * l^2; %antar masse [1 kg] på sattelitt
m_T = 5.972 * 10^24;

parameters = [M;
              m_T];


% "parameters" allows you to pass some parameters to the "SatelliteDynamics" function

time_final = 10; %Final time

% Simulate satellite dynamics
[time,statetraj] = ode45(@(t,x)SatelliteDynamics(t, x, parameters),[0,time_final],state);

% Here below is a template for a real-time animation
tic; % resets Matlab clock
time_display = 0; % time displayed

ScaleFrame = 5;   % Scaling factor for adjusting the frame size (cosmetic)
FS         = 15;  % Fontsize for text
SW         = 0.035; % Arrows size
while time_display < time(end)
    time_animate = toc; % get the current clock time
    % Interpolate the simulation at the current clock time
    state_animate = interp1(time,statetraj,time_animate);
     p = state_animate(1:3)';
    p = p*10^(-7);
    R = reshape(state_animate(4:12),3,3);
    omega = state_animate(16:18)';

    figure(1);clf;hold on
    

    % Use the example from "Satellite3DExample.m" to display your satellite
    MakeFrame(zeros(3,1),eye(3),ScaleFrame,FS,SW,'a', 'color', 'k')
    MakeFrame(p,R,ScaleFrame,FS,SW,'b', 'color', 'r')
    MakeArrow(p,R*omega,FS,SW,'$$\omega$$', 'color', [0,0.5,0])
    DrawRectangle(p,R ,'color',[0.5,0.5,0.5]);
    FormatPicture([0;0;2],0.5*[73.8380   21.0967   30.1493])

    if time_display == 0
        display('Hit a key to start animation')
        pause
        tic
    end
    time_display = toc; % get the current clock time
end
