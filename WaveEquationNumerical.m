%% Compute solution
% Computes the numerical solution to the wave equation.
% Supports both Dirichlet and Neumann boundary conditions.
% Also supports arbitrary initial conditions, and implements
% optional forcing wave motion in one corner of the domain,
% and an optional energy friction term.

clc, clear, close
L_x = 10;
L_y = 10;
T = 10;
N_x = 200;
N_y = 200;
N_t = 1000;
delta_x = L_x / N_x;
delta_y = L_y / N_y;
delta_t = T / N_t;

% Determines the energy loss due to friction. Must be >= 0.
mu = 0; 

c = 2;
stability_constant = (c*delta_t/delta_x)^2 + (c*delta_t/delta_y)^2;
disp("Stability constant: " +num2str(stability_constant))
sigma_x = c^2 * delta_t^2 / delta_x^2;
sigma_y = c^2 * delta_t^2 / delta_y^2;
x = (0:delta_x:L_x)';
y = (0:delta_y:L_y)';
u = zeros(N_x+1, N_y+1, N_t+1);

% Initial condition u(x,y,0) = g(x,y).
u(:,:,1) = g(x,y);

% Initial condition d_t u(x,y,0) = v(x). Does not affect the boundary.
V = v(x,y,L_x, L_y);
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


% Dirichlet boundary conditions:
alpha_1 = 0;
alpha_2 = 0;
alpha_3 = 0;
alpha_4 = 0;

u(:,1,:) = alpha_1;
u(:,end,:) = alpha_2;
u(1,:,:) = alpha_3;
u(end,:,:) = alpha_4;

% Computing the numerical solution
for k = 2:N_t
    for j = 2:N_y
        for i = 2:N_x
%           Neumann boundary counditions:
            f_1 = 0;
            f_2 = 0;
            f_3 = 0;
            f_4 = 0;
            
            u(i,1,k) = u(i,2,k) + + delta_t * f_1;
            u(i,end,k) = u(i,end-1,k) + delta_t * f_2;
            u(1,j,k) = u(2,j,k) + delta_t * f_3;
            u(end,j,k) = u(end-1,j,k) + delta_t * f_4;
            
            % Does does not change the solution, only makes it look better.
            u(1,1,k) = u(2,2,k);
            u(end,end,k) = u(end-1,end-1,k);
            u(end,1,k) = u(end-1,2,k);
            u(1,end,k) = u(1,end-1,k);
            
            % Generate wave:
            if (delta_t * k < T / 4)
                u(2,2,k) = sin(2 * 2 * pi * (delta_t * k) / (T / 4));
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

%% Plot the solution.
close
sleep_time = 0.1;

% Decreases the number of points to plot, increasing FPS.
total_time_points = 250;
t_skip = round(N_t / total_time_points);

total_x_points = 100;
total_y_points = 100;
x_skip = max(round(N_x / total_x_points),1); 
y_skip = max(round(N_y / total_y_points),1);

x_axis = (0:delta_x:L_x);
y_axis = (0:delta_y:L_y);

zmax = max(max(max(u)))+0.001;
zmin = min(min(min(u)))-0.001;
plotx = x_axis(1:x_skip:end);
ploty = y_axis(1:y_skip:end);
mySurf = surf(zeros(2,2),'FaceColor', 'interp', 'EdgeColor','black', 'EdgeAlpha', 0.5);

min_side_length = min(L_x, L_y);
pbaspect([L_x / min_side_length L_y / min_side_length 1])


axis ([0 L_x 0 L_y zmin zmax])
set(mySurf, 'linestyle', '-')
xlabel('x')
ylabel('y')
for i = 1:t_skip:N_t
    % Transpose u so that the x- and y-axes are shown correctly.
    plotz = u(1:x_skip:end,1:y_skip:end,i)';
    set(mySurf,'XData',plotx,'YData',ploty,'ZData',plotz, 'CData', plotz);
    drawnow;
    pause(delta_t*t_skip);
end


%% Functions
% Function determining initial velocity.
function v = v(x,y,L_x, L_y)
    v = zeros(length(x), length(y));
    for i = 1:length(x)
         for j = 1:length(y)
%             v(i,j) =  pulse(x(i)-L_x/2).*pulse(y(j)-L_y/2);
%             v(i,j) = 0;
        end
    end 
end

function g = g(x,y)
% Function determining initial position.
g = zeros(length(x), length(y));
end

function y = pulse(x)
sigma=5;
y = 1/(sqrt(pi)) * exp(-sigma*(x).^2);
end