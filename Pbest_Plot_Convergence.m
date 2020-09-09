%% Particle Swarm Optimization (PSO) on Levi Function (Minimization Problem)

clc;
clear all;
close all;

P = 100; % Total Number of Particles
max_iter = 100; 
% Since the range of x and y is [-10,10]
a = 0;
b = 20;
positions = (b-a).*rand(P,2)-10;
velocities = rand(P,2).*0;
particles = [positions,velocities];

% In the 'particles' array, the 1st two columns represent 'positions' and
% the next two represent 'velocities'

for i=1:P 
    % Deciphering the fitness value (kept in the 5th column)
    particles(i,5)=sin(3*pi*particles(i,1))^2 + (particles(i,1)-1)^2 * (1 + (sin(3*pi*particles(i,2)))^2) + (particles(i,2)-1)^2 * (1 + (sin(2*pi*particles(i,2)))^2); 
end

% Sorting the rows with respect to the column consisting the fitness value
particles = sortrows(particles,5);
 
Pbest = [particles(:,1:2),particles(:,5)]; % Initially Pbest(i) is the copy of positions(i)
Gbest = [particles(1,1:2),particles(1,5)]; % This is the best fit particle in the 1st iteration
 
Gbest_list = Gbest;
new_variable = Pbest;

% keeping the Pbest columns also in the main 'particles' array (6th, 7th & 8th columns)
particles = [particles, Pbest];

colNames = {'x','y','v_x','v_y','f(x,y)','Pbest_x','Pbest_y','Pbest_func'};
sTable = array2table(particles,'VariableNames',colNames);

% Reference:
% https://www.sciencedirect.com/topics/computer-science/particle-swarm-optimization :
% used for choosing these coefficients
alpha = rand(1); % Inertia Coefficient
c1 = rand(1); % Personal Acceleration Coefficient
c2 = rand(1); % Social Acceleration Coefficient

%alpha = 0.6329;
%c1 = 0.7927;
%c2 = 0.6930;

%alpha = 1.2;
%c1 = 1;
%c2 = 2;

% If these three coefficients are within the range [0,1], the solutions
% obtained very fastly. 
 
% Iterating it for 100 times
for j=1:max_iter
    for k=1:P
       r1 = rand(1);
       r2 = rand(1);
       
       % Updating the velocities
       % Velocity (x)
       particles(k,3) = alpha * particles(k,3) + c1 * r1 * (particles(k,6)-particles(k,1)) + c2 * r2 * (Gbest(1,1)-particles(k,1));
       % Velocity (y)
       particles(k,4) = alpha * particles(k,4) + c1 * r1 * (particles(k,7)-particles(k,2)) + c2 * r2 * (Gbest(1,2)-particles(k,2));
       
       %Updating the positions
       % Position x
       particles(k,1) = particles(k,1) + particles(k,3);
       % Position y
       particles(k,2) = particles(k,2) + particles(k,4);
       
       %Ensuring the range of 'x' and 'y' within [-10,10] range
       if (particles(k,1)<-10)
           particles(k,1)=-10;
       elseif (particles(k,1)>10)
           particles(k,1)=10;
       end
           
       if (particles(k,2)<-10)
           particles(k,2)=-10;
       elseif (particles(k,2)>10)
           particles(k,2)=10;
       end
       
       % Evaluating the fitness value
       particles(k,5)=sin(3*pi*particles(k,1))^2 + (particles(k,1)-1)^2 * (1 + (sin(3*pi*particles(k,2)))^2) + (particles(k,2)-1)^2 * (1 + (sin(2*pi*particles(k,2)))^2);
    end
    
    % Sorting the rows w.r.t the fitness value
    particles = sortrows(particles,5);
    
    % Updating the Global best 
    if ((Gbest(1,3))>particles(1,5))
           Gbest(1,1:2)=particles(1,1:2);
           Gbest(1,3)=particles(1,5);
    end
    
    % Appending the Gbest values every iteration
    Gbest_list=[Gbest_list;Gbest];
    
    % Updating the Pbest
    % Since this is a minimization problem, if the Pbest is greater than 
    % the current fitness value of a particle, then the Pbest value for 
    % that particle is updated
    for l=1:P
       if (particles(l,8)>particles(l,5))
           particles(l,6:7)=particles(l,1:2);
           particles(l,8)=particles(l,5);
       end
    end
    new_variable = [new_variable;particles(:,6),particles(:,7),particles(:,8)]; % for plotting convergence of Pbest
end

Answer = Gbest_list(max_iter+1,:)

sTable = array2table(particles,'VariableNames',colNames);


%% Plotting the convergence of Pbest

[X,Y] = meshgrid(-10:0.5:10);
Z = sin(3*pi.*X).^2 + ((X-1).^2).*(1+(sin(3*pi.*Y).^(2))) + ((Y-1).^2).*(1+((sin(2*pi.*Y)).^(2)));
s = surf(X,Y,Z,'FaceAlpha',0.1);
zlim([-5 200]);

hold on;

set(gca,'XLim',[-10 10],'YLim',[-10 10]);
view(43,24);

title('Pbest - 1st Iteration');
xlabel('X');
ylabel('Y');
zlabel('Levi Function');
pause(0.01);

for i=1:100
   a = 'Iteration ' + string(i);
   title(a);
   num = (i-1)*100;
   head = scatter3(new_variable(num+1:num+100,1),new_variable(num+1:num+100,2),new_variable(num+1:num+100,3),'filled','b');
   tail = scatter3(Gbest_list(i,1),Gbest_list(i,2),Gbest_list(i,3),'filled','r')
   drawnow
   F(i) = getframe(gcf);
   pause(0.1);
   delete(head)
   delete(tail)
end

video = VideoWriter('yobaby.avi','Uncompressed AVI');
video.FrameRate = 1;
open(video)
writeVideo(video,F);
close(video)

