function [min_dist, best_route_per_robot, best_dist_travelled_robot, time_elapsed, iter_tests] = MTSP_RANDOM(PoI,dist_matrix,current_robots,robots,clr,PoILon,PoILat,max_distance,weight_F2,flag_for_checking_time,max_time,iter_tests,indiv_results,indiv_hist)
%% Random function - calculates and evaluates a possible PoI and robot combinations per iteration

fprintf('Begin Random! Setting up... \n');

%% Setting up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

best_dist_travelled_robot = 1e100*ones(1,current_robots);
min_dist = inf; %minimal distance
acum_dc = []; %to store the minimal distances
time_to_iter = []; %vector to save the iteration of the best solution
x_ticks = 1;
results = [];
flag_for_reached_time_limit = 0;
j = 1;

fprintf('Computing for %d robots! \n', current_robots)
tic;

%% Loop starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for max_iteration = iter_tests % number of iterations for the loop, calculate minimal distance
    if flag_for_reached_time_limit == 1
		break;
    end
    iteration = 0;
    while iteration < max_iteration
        seq_robots = ones(1,PoI-1); %actual sequence of robots to visit the panel
		total_elapsed_time = toc;
        if flag_for_checking_time %loop for time calculation
            if iteration > x_ticks(end)
                time_to_iter(end+1) = iteration;
                x_ticks(end+1) = x_ticks(end)*10;
                results(end+1) = min_dist;
            end
            if total_elapsed_time > max_time
			    flag_for_reached_time_limit = 1;
			    break;
            end
        end

        iteration = iteration + 1;
        route = randperm(PoI-1)+1; %taking a random route to explore

        %calculating the distance considering random combinations
        dist_travelled_robot = zeros(1,current_robots); %distance for each robot
        total_current_distance = 0; %total  current distance
        current_robot_position = ones(1,current_robots); %last panel visited by each robot
        prox = 0;  %flag variable

        % randomly change the robots in seq_robots until
        % the amount of different robots is greater that current_robots
        while length(unique(seq_robots)) < current_robots
            for i=1:length(seq_robots)
                seq_robots(i) = randi(current_robots);
            end
        end

        for k = 1:(PoI-1)
            %accumulating distance
            total_current_distance = total_current_distance + dist_matrix(current_robot_position(seq_robots(k)),route(k)); 

            %distance for current robot
            dr = dist_matrix(current_robot_position(seq_robots(k)),route(k));
            dist_travelled_robot(seq_robots(k)) = dist_travelled_robot(seq_robots(k)) + dr;

            %robot distance is bigger then max allowed distance for
            %robot? go to the next robot combination
            if (dist_travelled_robot(seq_robots(k)) + dist_matrix(1,route(k))) > max_distance
                prox = 1;
                break
            end
            current_robot_position(seq_robots(k)) = route(k);
        end

        % case the distances for each robot were smaller then de max distance allowed,
        % evaluate if the distance is the mininal
        if prox == 0
            for k = 1:current_robots
                total_current_distance = total_current_distance + dist_matrix(1,current_robot_position(k)); %SOMA DISTANCIA DO 1 Á ULTIMA CASA PERCORRIDA    
                dist_travelled_robot(k) = dist_travelled_robot(k) + dist_matrix(1,current_robot_position(k)); 
            end

            acum_dc =[acum_dc total_current_distance];  %storing the distance


			max_dist_robot = max(dist_travelled_robot);
			min_dist_robot =  min(nonzeros(dist_travelled_robot));
			best_second_term = 1e100;
			if ~isequal(1e100*ones(1,current_robots), best_dist_travelled_robot)
				best_second_term = weight_F2*(max(best_dist_travelled_robot)...
					- min(nonzeros(best_dist_travelled_robot)));
			end

            if (1-weight_F2)*total_current_distance + weight_F2*(max_dist_robot - min_dist_robot)... %fitness function comparison
			< (1-weight_F2)*min_dist  + best_second_term
				best_dist_travelled_robot = dist_travelled_robot;
                min_dist_robot = dist_travelled_robot;
                min_dist = total_current_distance;    % Minimum total route distance
                best_route = route; % Vector that contains the order of panels inside the best route
                robot_routes = seq_robots; %storing the best robot sequence visit
            end
        end
    end 
    results(end+1) = min_dist;
end
time_elapsed = toc;

fprintf('Finished! \nPlotting Graph \n')

%% Creates vector with best route solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

best_route_per_robot = ones(current_robots,PoI+1);
for k = 1:current_robots      
    %route for a robot
    i=1;
    for x = 1:(PoI-1)
        if robot_routes(x) == k
            best_route_per_robot(k,i) = best_route(x);
            i=i+1;
        end  
    end
end

%% Plot every best result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if indiv_results == 1
    figure()
    PoIcoords = [PoILon;PoILat];
    plot_results_no_charge(PoIcoords',change_routes_format(best_route_per_robot),clr);
    title(sprintf('Total Distance = %1.1f, Iteration = %d \n',min_dist,iteration))
    daspect([1 1 1])
end

%% Plots history of results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if indiv_hist == 1
    figure_id = 100 + current_robots;
    figure(figure_id)
    sgtitle('MTSP Random results history: total dist /  number of iteration ')
    
    if flag_for_checking_time
        semilogx(x_ticks,results,'DisplayName',"" + robots + " robots","LineWidth",3,'Marker','*','Color', clr(k,:));
	    xlim([1 max(time_to_iter)])
        xlabel("Number of iterations [#]",'FontSize',12)
    else
	    semilogx(iter_tests,results,"LineWidth",2,'Marker','*','Color', clr(k,:));
	    xlim([1 (max(iter_tests))])
        xlabel("Number of iterations [#]",'FontSize',12)
    end
    legend("" + robots + " robots")
    fprintf("")
    ylabel("Total dist [m]")
    grid on
    hold on
end

figure(105)
subplot(2,ceil(robots/2),current_robots);
set(gca,'XScale','log')
hold on
if flag_for_checking_time
    semilogx(x_ticks,results,'DisplayName',"" + current_robots + " robots","LineWidth",3,'Marker','*','Color', clr(k,:));
	xlim([1 max(time_to_iter)])
    xlabel("Number of iterations [#]",'FontSize',12)
else
	semilogx(iter_tests,results,"LineWidth",2,'Marker','*','Color', clr(k,:));
	xlim([1 (max(iter_tests))])
    xlabel("Number of iterations [#]",'FontSize',12)
end
legend("" + current_robots + " robots")
fprintf("")
ylabel("Total dist [m]")
grid on
hold on
end