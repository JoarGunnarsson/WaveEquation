clear, clc
% The analytic equation to the wave equation with Dirichlet boundary
% conditions, with a constant initial velocity across the entire domain,
% in a square domain.

T = 10;
N_time = 250;
delta_t = T/N_time;

a = 1;
c = 1;
K = 5;
N_space = 100;
delta_x = a / N_space;
delta_y = a / N_space;

iter = 10; % The number of terms in the sum.

u = zeros(N_space+1, N_space+1, N_time);
for k = 1:N_time
    for x = 1:N_space+1
        for y = 1:N_space+1
            for i = 1:iter
                for j = 1:iter
                    u(x,y,k) = u(x,y,k) - K * (1 - (-1)^i) / i * (1 - (-1)^j) / j * 1/(pi^2*i^2/a^2 + pi^2*j^2/a^2)^0.5 * sin(c * (pi^2*i^2/a^2 + pi^2*j^2/a^2)^0.5 * k * delta_t) * sin(pi*i*x*delta_x / a) * sin(pi*j*y*delta_y / a);
                end
            end
        end
    end
end

disp("Done calculating.")
%% Plotting soluton

close

% Decreases the number of points to plot, increasing FPS.
x_skip = max(round(N_space / 100),1);
y_skip = max(round(N_space / 100),1);
x_axis = (0:delta_x:a);
y_axis = (0:delta_y:a);

zmax = max(max(max(u)));
zmin = min(min(min(u)));
time_wait_factor = 1;
plotx = x_axis(1:x_skip:end);
ploty = y_axis(1:y_skip:end);
plotz = u(1:y_skip:end,1:x_skip:end,1);
mySurf = surf(zeros(2,2),'FaceColor', 'red','EdgeColor','black');
axis ([0 a 0 a zmin zmax])
set(mySurf, 'linestyle', '-')
for i = 1:N_time
    plotz = u(1:y_skip:end,1:x_skip:end,i);
    set(mySurf,'XData',plotx,'YData',ploty,'ZData',plotz);
    
    
    pause(delta_t * time_wait_factor)
end



