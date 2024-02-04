%% Compute solution
clc, clear, close
L_x = 20;
L_y = 10;
T = 10;
N_x = 50;
N_y = 50;
N_t = 200;
delta_x = L_x / N_x;
delta_y = L_y / N_y;
delta_t = T / N_t;
mu = 0; % Determines the energy loss due to friction. Must be >= 0.

c = 2;
stability_constant = (c*delta_t/delta_x)^2 + (c*delta_t/delta_y)^2;
disp("Stability constant: " +num2str(stability_constant))
sigma_x = c^2 * delta_t^2 / delta_x^2;
sigma_y = c^2 * delta_t^2 / delta_y^2;
x = (0:delta_x:L_x)';
y = (0:delta_y:L_y)';
u = zeros(N_x+1, N_y+1, N_t+1);

% Initial condition u(x,y,0) = 0.
u(:,:,1) = 0;

% Initial condition d_t u(x,y,0) = v(x). Does not affect the boundary
V = v_0(x,y,L_x, L_y);
ghost_cells = zeros(N_x+1, N_y+1);


for j = 2:N_y
    for i = 2:N_x
        t_part = 2 * u(i,j,1) + V(i,j) * 2 * delta_t;
        x_part = sigma_x * (u(i+1,j,1) - 2*u(i,j,1) + u(i-1,j,1));
        y_part = sigma_y * (u(i,j+1,1) - 2*u(i,j,1) + u(i,j-1,1));
        u(i,j,2)  = t_part + x_part + y_part;
        u(i,j,2) = u(i,j,2) - mu * delta_t / 2 * V(i,j) * 2 * delta_t / (1 + mu * delta_t / 2);
        u(i,j,2) = u(i,j,2) / ( (2 - mu * delta_t / 2) / (1 + mu * delta_t / 2));
    end
end


% Dirichlet boundary conditions.
u(1,:,:) = 0;
u(end,:,:) = 0;
u(:,1,:) = 0;
u(:,end,:) = 0;

% Computing the numerical solution
for k = 2:N_t
    for j = 2:N_y
        for i = 2:N_x
%           Neumann boundary counditions dx u and/or dy u = 0:            
            u(1,j,k) = u(2,j,k);
            u(end,j,k) = u(end-1,j,k);
            u(i,1,k) = u(i,2,k);
            u(i,end,k) = u(i,end-1,k);
            
            % Does does not change the solution, only makes it look better.
            u(1,1,k) = u(2,2,k);
            u(end,end,k) = u(end-1,end-1,k);
            u(end,1,k) = u(end-1,2,k);
            u(1,end,k) = u(1,end-1,k);
            
            % Generate wave
             if (delta_t * k < T / 4)
                 u(N_x / 2,N_y / 2,k) = sin(2 * 2 * pi * (delta_t * k) / (T / 4));
             end
            
            t_part = 2 * u(i,j,k) - u(i,j,k-1);
            x_part = sigma_x * (u(i+1,j,k) - 2*u(i,j,k) + u(i-1,j,k));
            y_part = sigma_y * (u(i,j+1,k) - 2*u(i,j,k) + u(i,j-1,k));
            u(i,j,k+1) = t_part + x_part + y_part;
%           Introduces an approximative friction loss term:
            u(i,j,k+1) = (u(i,j,k+1) + mu * u(i,j,k-1) * 1/2 * delta_t) / (1 + mu * delta_t / 2);
        end
    end
end

%% Plot any solution.
close
sleep_time = 0.1;

total_time_points = 250;
t_skip = round(N_t / total_time_points);

x_skip = max(round(N_x / 1000),1); %Gör så att man inte plottar hela u, ökar prestanda.
y_skip = max(round(N_y / 1000),1); %Gör så att man inte plottar hela u, ökar prestanda.

x_axis = (0:delta_x:L_x);
y_axis = (0:delta_y:L_y);

zmax = max(max(max(u)))+0.001;
zmin = min(min(min(u)))-0.001;
plotx = x_axis(1:x_skip:end);
ploty = y_axis(1:y_skip:end);
mySurf = surf(zeros(2,2),'FaceColor', 'interp', 'EdgeColor','black');

min_side_length = min(L_x, L_y);
pbaspect([L_x / min_side_length L_y / min_side_length 1])
axis ([0 L_x 0 L_y zmin zmax])
set(mySurf, 'linestyle', '-')
xlabel('x')
ylabel('y')
for i = 1:t_skip:N_t
    % Transpose u so that the x- and y-axes are shown correctly.
    plotz = u(1:x_skip:end,1:y_skip:end,i)';
    set(mySurf,'XData',plotx,'YData',ploty,'ZData',plotz);
    
    pause(delta_t)
end


%% Functions
% Function determining initial velocity of the domain.
function result = v_0(x,y,L_x, L_y)
    result = zeros(length(x), length(y));
    for i = 1:length(x)
         for j = 1:length(y)
%             result(i,j) = 0.1*sin(x(i)* 2 * pi / L_x) * sin(y(j) * 2 * pi / L_y);
%             result(i,j) = result(i,j) +pulse(x(i)-L_x/6).*pulse(y(j)-L_y/6);
%             result(i,j) = result(i,j) + pulse(x(i)-L_x/6).*pulse(y(j)-5*L_y/6);
%             result(i,j) = result(i,j) + pulse(x(i)-5*L_x/6).*pulse(y(j)-L_y/6);
%             result(i,j) = result(i,j) + pulse(x(i)-5*L_x/6).*pulse(y(j)-5*L_y/6);
%             result(i,j) = result(i,j) + pulse(x(i)-L*0.9).*pulse(y(j)-L_y/2);
%             result(i,j) = 1*pulse(x(i)-L_x/2).*pulse(y(j)-L_y/2);
             result(i,j) = 0;
        end
    end 
end