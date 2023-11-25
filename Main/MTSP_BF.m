function [min_dist, best_route_per_robot_cases,time_elapsed, case_dc] = MTSP_BF(PoI,dist_matrix,robots,clr,PoILon,PoILat,max_distance,weight_wdm,disp_time_bf_milp)
% Brute force function - calculates and evaluates all possible PoI and robot combinations

fprintf('Begin Brute Force! Setting up... \n');

%% Setting up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PoIcoords = [PoILon;PoILat];
exp_routes = perms(2:PoI);  	%Finds all possible PoI combinations
comb_robot = distribution(robots, PoI-1); % uses function to find all the robot combinations
robot_routes_cases = zeros(size(comb_robot));
total_robot_combinations = size(comb_robot,1);
min_dist = inf; 	%minimal distance
min_dist_robot = zeros(1,robots); %minimal distance for each robot in the best solution 
acum_dc = []; 	%to store the minimal distances
case_dc = zeros(1,robots); %to store the minimal distances but per each number of robot available
case_dc(:) = inf;
best_dist_travelled_robot = zeros(1,robots); %to store the distances for each robot for the current best case 
time_elapsed = zeros(1,robots); % vector to store the time it took to calculate said minimal distances per each number of robot available
seq_route = ones(1,PoI-1);
best_route_per_robot_cases = ones(robots,PoI+1,robots);
best_route_cases = ones(robots,PoI-1);
fprintf('Setup done! Computing... \n')

%% Loop starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:size(exp_routes,1) % loop to calculate minimal distance

    if mod(j, 100) == 0
        j
    end

    zeros_count1 = (robots-1);    
    route = exp_routes(j,:); %taking a route to explore
    index_dc = 1;

    for m = 1:total_robot_combinations %calculating the distance considering each combination of robots

        tic;
        dist_travelled_robot = zeros(1,robots); %distance for each robot
        dc = 0;  %total current distance
        prox = 0;  %flag variable 
        ind = 1;
        zeros_count2 = 0;
        
        for r = 1:robots
    
            f = comb_robot(m,r); %select robot combination
    
            if comb_robot(m,r) == 0 
    
                zeros_count2 = zeros_count2 + 1; 
    
            end
    
            if f > 0
    
                [seq_route(ind:ind+f-1),dr]= tsp(route(ind:ind+f-1),dist_matrix);                
                ind = ind+f;
                dist_travelled_robot(r) = dr;  %distance for current robot 
                dc = dc + dr; %accumulating distance   
    
                %robot distance is bigger then max allowed distance for
                %robot? go to the next robot combination
                if (dr > max_distance)
                    prox = 1;
                    break
                end
            end
        end
    
        if  zeros_count2 < zeros_count1
            time_elapsed(index_dc) = toc + time_elapsed(index_dc);
            tic;
            index_dc = index_dc + 1;
            zeros_count1 = zeros_count2;
        end
    
        % case the distances for each robot were smaller then de max distance allowed,
        % evaluate if the distance is the minimal
        if prox == 0
    
            max_dist_robot = max(dist_travelled_robot);
            min_dist_robot = min(nonzeros(dist_travelled_robot));
            best_second_term = 0;
            if ~isequal(zeros(1,robots), best_dist_travelled_robot)
                best_second_term = weight_wdm*(max(best_dist_travelled_robot) - min(best_dist_travelled_robot));
            end
            if (1-weight_wdm)*dc + weight_wdm*(max_dist_robot - min_dist_robot)... % fitness function comparison
                    < (1-weight_wdm)*case_dc(index_dc)  + best_second_term
    
                best_dist_travelled_robot = dist_travelled_robot;
                best_route_cases(index_dc,:) = seq_route;
                best_route_per_robot_cases(:,:,index_dc) = ones(robots,PoI+1);
                robot_routes_cases(index_dc,:) = comb_robot(m,:);
                case_dc(index_dc) = dc;	
                i = 1;
    
                for r = 1:robots
                    for x = 1:robot_routes_cases(index_dc,r)
    
                        best_route_per_robot_cases(r,x,index_dc) = seq_route(i);            		
                        i=i+1;
    
                    end
                end
            end
    
            acum_dc = [acum_dc dc];
    
            if min_dist > dc               
    
                min_dist = dc;   % Minimum total trip distance
                best_trip = seq_route; % Vector that contains the order of panels inside the best trip
                robot_trips = comb_robot(m,:);  %storing the best robot sequence visit
    
            end
        end
    end

    time_elapsed(robots) = toc + time_elapsed(robots);

end 

fprintf('Finished! \nPlotting Graph... \n')

%% Plot with best trip solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if disp_time_bf_milp == 0
    figure
    for r = 1:robots
        subplot(2,ceil(robots/2),r);
        plot_results_no_charge(PoIcoords',change_routes_format(best_route_per_robot_cases(:,:,r)),clr);
        daspect([1 1 1])
        title(sprintf('Total Distance = %1.1f',case_dc(r)),"FontWeight","normal");
    end
    sgtitle('MTSP BF')
    display(case_dc);
    pad=10;
    if(~isempty(acum_dc))
        figure
        bar(acum_dc)
        xlim([0 length(acum_dc)])
        ylim([(min(acum_dc)-pad) (max(acum_dc)+pad)])
        ylabel('Total dist [m]')
        xlabel('Combination [#]')
        title('MTSP BF Results: total dist / number of combinations',"FontWeight","normal")
        hold on
    end
end

%% Extra functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dist = distribution(robots, panels) %function to calculate the robot combinations, gives number of panels for each robots
    
        max_value = panels;
        total = max_value^robots;
        comb = uint8(zeros(1, robots));   
        dist = [];
        
        for i=1:total                   
            for j = robots:-1:1
                if (sum(comb)==panels) && (sum(comb(1:end-1)<= comb(2:end))== (robots-1))
                    dist = [dist; comb];
                end
                comb(j) = comb(j) + 1;
                
                if comb(j) > max_value 
                    comb(j) = 0;
                else
                    break
                end
            end
        end    
    end

    function [seq, min] = tsp(panels, dist_matrix) %function to calculate the best distance for TSP
        if panels<2
            seq = panels;
            min = 2*dist_matrix(1,panels);
        else
            trips = perms(panels);
            min = inf;
            indx = -1;
            dcs =  zeros(1,size(trips,1));
            
            for z=1: size(trips,1)
                %saving the distance from 1 to the first PoI
                current_dist = dist_matrix(1, trips(z,1));
                %acummulation the distance
            
                for h=2:length(panels)
                    current_dist = current_dist + dist_matrix(trips(z,h-1), trips(z,h));
               
                    %avoiding trips  bigger then the minimal
                    if current_dist>min
                        break
                    end
                    
                end
                
                %adding the return distance to depot      
                current_dist = current_dist + dist_matrix(trips(z,end), 1);            
                dcs(z) = current_dist;
                if current_dist< min
                    min = current_dist;
                    indx = z; 
                end                        
            end                    
            seq = trips(indx,:);        
        end
    end
end