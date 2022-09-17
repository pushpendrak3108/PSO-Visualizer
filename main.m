% Purpose: Use the particle swarm optimization PSO algorithm to find the maximum value of the three-dimensional function
clc; clear; close all;

num = 128;% number of individuals in the flock
c1 = 2; c2 = 2; %c1 c2 = learning factor (non-negative constant)
w = 0.8; % inertia factor (non-negative)
alpha = 0.4; % constraint factor = weight to control speed
v_limit = 0.2;% speed limit

%Step2=Initialize the solving function
[fun_x, fun_y, fun_z] = peaks; % PEAKS is a function of two variables, obtained by translating and scaling Gaussian distributions, which is useful for demonstrating MESH, SURF, PCOLOR, CONTOUR, etc.
figure(1);
hold on;
meshc(fun_x, fun_y, fun_z);

%Step3=Use PSO to solve the maximum value of the function
particle = zeros(num, 5);% a group of particles per row = [number position x position y velocity x velocity y]
pi = zeros(num, 3); % individual historical optimal each row = [number position x position y] (each individual has a historical best and historical best with memory = best so far, including the past status)
pg = zeros(1, 2); % group history is optimal each row = [position x position y] (there is only one for the entire particle swarm)
pg(1) = -3 + 6*rand(1,1); %x=-3~3
pg(2) = -3 + 6*rand(1,1);%y=-3~3
% Give the initial value of the particle swarm (random position and speed of each particle)
for i = 1: 1: num
    % per particle
    particle(i,1) = i; %Number
    particle(i,2) = -3+6*rand(1,1); %x=-3~3
    particle(i,3) = -3+6*rand(1,1); %y=-3~3
    particle(i,4) = -v_limit+v_limit*2*rand(1,1);%vx=-v_limit~v_limit
    particle(i,5) = -v_limit+v_limit*2*rand(1,1);%vy=-v_limit~v_limit
    % Individual History Best
    pi(i,1) = particle(i,1); %Number
    pi(i,2) = particle(i,2);%x
    pi(i,3) = particle(i,3);%y
    % group historical optimal initial value
    if GetFitness( particle(i,2:3) ) > GetFitness( pg )
        pg = particle(i,2:3);
    end
end
% Display of the initial situation
figure(1);
hold on;
plot_x = particle(:,2);
plot_y = particle(:,3);
plot_z = FuncCalculate(plot_x, plot_y);
plot3(plot_x , plot_y , plot_z , '*r');% particle position
axis([-3,3 , -3,3]);
hold off;

% prepare avi animation
aviObj = VideoWriter('Particle Swarm Optimizer.avi');% save as avi
aviObj.FrameRate = 5;
open(aviObj);

for cycle=1:1:100% particle swarm movement times
    tmp_pg = pg; % group optimal
% The movement of each particle should be parallel, so the refreshed group history optimization should not affect the movement of this round
% For each particle PSO motion is based on pg
% Refresh pg after all particle motions are over

    for i = 1: 1: num% each particle
        %1) Particle operation (see: Particle Swarm Optimization Algorithm Li Aiguo, Qin Zheng, Bao Fumin, He Shengping Equations 3 and 4)
        v_id = particle(i,4:5); %i particle velocity=(v_x v_y)
        x_id = particle(i,2:3); %i particle position=(x_x x_y)
        p_id = pi(i, 2:3); %the individual historical optimal position of the i particle
        
        r1 = rand(1,1);%r1,r2 are random numbers between [0,1]
        r2 = rand(1,1);
        
        v_id = w*v_id + c1*r1*(p_id-x_id) + c2*r2*(pg-x_id);
        x_id = x_id + alpha*v_id;
        
        %2) Adjustment for out-of-bounds particles: on the reflection boundary = the speed does not change the direction is reversed
        if ( x_id(1)>3 && v_id(1)>0 ) || ( x_id(1)<-3 && v_id(1)<0 )% x+ or x- direction out of bounds
            x_id = x_id - alpha*v_id;% restore the previous position
            v_id(1) = -v_id(1); %Invert the speed in the direction of x
            x_id = x_id + alpha*v_id;% and then move
        end
        if ( x_id(2)>3 && v_id(2)>0 ) || ( x_id(2)<-3 && v_id(1)<0 )%y+ or y-direction out of bounds
            x_id = x_id - alpha*v_id;% restore the previous position
            v_id(2) = -v_id(2); %Inverse the speed in the y direction
            x_id = x_id + alpha*v_id;% and then move
        end
        
        %3) Particle motion % refresh the position and speed of the particle
        particle(i,2:3) = x_id;
        particle(i,4:5) = v_id;
        % refresh individual history optimal
        if GetFitness( particle(i,2:3) ) > GetFitness( pi(i,2:3) )
            pi(i,2:3) = particle(i,2:3);
        end
        % refresh group history best
        if GetFitness( particle(i,2:3) ) > GetFitness( tmp_pg )
            tmp_pg = particle(i,2:3);%tmp_pg does not affect the subsequent particle judgment in this round
        end
    end
    pg=tmp_pg;
    
% Visualize the current result
    figure(2);
	%Function image
    meshc(fun_x, fun_y, fun_z);
    hold on;
% distribution of each particle
    plot_x = particle(:,2);
    plot_y = particle(:,3);
    plot_z = FuncCalculate(plot_x, plot_y);
    plot3( plot_x , plot_y , plot_z , '*r');% particle position
    axis([-3,3 , -3,3]);
    hold off;
    % write avi animation
    writeVideo(aviObj,getframe(gcf));
end
% save avi animation
close(aviObj);

% The final result is displayed
figure(3);
meshc(fun_x, fun_y, fun_z);% function image
hold on;
plot_x=particle(:,2);
plot_y=particle(:,3);
plot_z=FuncCalculate(plot_x, plot_y);
plot3( plot_x , plot_y , plot_z , '*r');% particle position
axis([-3,3 , -3,3]);
hold off;