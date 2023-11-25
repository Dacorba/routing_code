function varargout = MTSP_GAC_MSD(xy,dmat,charging_points,robots,min_route,max_distance,weight_wdm,pop_size,num_iter,use_complex,show_prog,show_res,clr)
%   MTSP_GAC_MSD, MTSP GA with fixed charging points, single/multiple depot
%   Finds a near best solution to a variation of the MTSP by setting up a 
%   GA to search for the shortest route (least distance needed for the robots 
%   to travel to each PoI exactly once and return to their starting locations). 
%   The robots start and end at the depot, which is a charging point. Although 
%   more charging points (which are not depots) can be included by the user in 
%   the Main script. This algorithm is based on Elad Kivelevitch's MDMTSPV_GA, 
%   but transforms the multiple depots into charging points, keeping a single
%   depot. Although, there is the option of multiple depot: in this case,
%   every charging point is also a depot.
% 
% Summary:
%     1. Each robot travels to a unique set of PoI and completes the
%        route by returning to the depot he started from
%     2. Each PoI is visited by exactly one robot
%     3. The total distance travelled by the robots used is minimized 
%        during the algorithm, while also minimizing the difference
%        between, the biggest distance travelled by the robots,
%        and the smallest distance travelled by the robots.
%     4. There are charging points which allow charging for the robots.
%     5. These charging points' location are defined in the Main script by 
%        the user.
%
% Note: The depot is taken to be the first xy point
%
% Input:
%     xy (float) is an Nx2 matrix of PoI locations, where N is the number
%               of PoI
%     dmat (float) is an NxN matrix of PoI-to-PoI distances or costs
%     charging_points (vector) extra locations for the robots to charge
%               mid-route
%     robots (scalar integer) is the number of robots to visit the PoI
%     min_route (scalar integer) is the minimum number of PoIs for each
%				robots, NOT including the depot
%     max_distance (scalar integer) is the maximum route length for each robot
%     weight_wdm (float) determines how much importance is given in the balance to 
%               minimize the distance between the highest and lowest values,
%               and minimizing the total distance
%     pop_size (scalar integer) is the size of the population (should be divisible by 8)
%     num_iter (scalar integer) is the number of desired iterations for the algorithm to run
%     use_complex (scalar boolen 0/1) is the flag wether to use complex mutation operators or not
%     show_prog (scalar logical) shows the GA progress if true
%     show_res (scalar logical) shows the GA results if true
%     clr (3 x robots float matrix) is the set of color to be used for
%               drawing the routes, where each line represent one robot
% Output:
%     opt_rte (integer array) is the best route found by the algorithm
%     min_dist (scalar float) is the total distance traveled by the robots
%     opt_iter (scalar int) is the number of iterations until the optimal solution has been found
%	  opt_time (scalar float) is the time in milliseconds until the optimal solution has been found
%     dist_history (float array) is the history of distances of best found solutions
%     cases_min (robots x 1 vector) is the minimum total distance for each
%               robot combination
%     cases_min_F (robots x 1 vector) is the minimum total fitness for each
%               robot combination
%     cases_routes_min (robots x 1 cell array) is the route that has the minimum total fitness for each
%               robot combination
%     cases_opt_time (robots x 1 vector) is the time it took to find the best
%               solution for each robot combination
%     cases_opt_iter(robots x 1 vector) is the iteration number when it
%               found the best solution
%     cases_each_robot_dist (robots x robots matrix) is the distance
%               travelled by each robot for each robot combination
%
% Route/Breakpoint Details:
%     The algorithm uses a data structure in which RTE lists the PoI and 
%     charging points in a route and BRKS lists break points that divide 
%     RTE between robots. This part of the algorith works exactly as 
%     described before on the MTSP_GA_SD, but note that if the robots had to 
%     charge mid-route, and two extra charging points were given, it could 
%     look like the following:
%         For robot 1: [0 5 -1 6 9 0] - from depot along the route and back
%         For robot 2: [0 1 4 2 0 8 0] - from depot along the route and back
%         For robot 3: [0 10 3 -1 7 0] - from depot along the route and back
%     That is, the negative numbers will represent the visit to the
%     aforementioned extra charging points.
% Author: Daniel Correia Batista
% Email: daniel.batista@tecnico.ulisboa.pt
% Based on: Elad Kivelevitch's MDMTSPV_GA (see MATLAB Central for download)
% Release: 1.0
% Release Date: 28/10/2023

fprintf('Begin MTSP w/ Charging Points \n');

%% Setting up%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the distance matrix for the charging points====================
charging_points = [xy(1,1:2); charging_points]; % isolate the charging points from the PoI, includes the depot
xy = xy(2:end,1:2);
cmat = zeros(size(xy,1),size(charging_points,1)); % matrix for the distances from charging points to PoI
for p = 1:size(xy,1)
    for h=1:size(charging_points,1)
        cmat(p,h) = norm(xy(p,:)-charging_points(h,:));
    end
end
%==========================================================================
merging_prob = 0.3;
%==========================================================================
% Verify Inputs============================================================
[N, ~] = size(xy);
[nr,nc] = size(dmat);
if N ~= nr || N ~= nc
    error('Invalid XY or DMAT inputs!')
end
n = N;
%==========================================================================
% Sanity Checks============================================================
robots = max(1,min(n,round(real(robots(1)))));
min_route = max(1,min(floor(n/robots),round(real(min_route(1)))));
pop_size = max(8,8*ceil(pop_size(1)/8));
num_iter = max(1,round(real(num_iter(1))));
show_prog = logical(show_prog(1));
show_res = logical(show_res(1));
%==========================================================================
% Initializations for Route Break Point Selection==========================
num_brks = robots-1;
dof = n - min_route*robots; % degrees of freedom
addto = ones(1,dof+1);
for k = 2:num_brks
	addto = cumsum(addto);
end
cum_prob = cumsum(addto)/sum(addto);
% =========================================================================
% Initialize the Populations===============================================
pop_rte = zeros(pop_size,n);          % population of routes
pop_brk = zeros(pop_size,num_brks);   % population of breaks
for k = 1:pop_size
	pop_rte(k,:) = randperm(n);
	pop_brk(k,:) = randbreaks();
end
% =========================================================================
% Run the GA===============================================================
global_min      = Inf;
cases_min		= Inf * ones(1,robots);
cases_min_F	= Inf * ones(1,robots);
cases_routes_min = cell(1,robots);
cases_opt_time = zeros(1,robots);
cases_opt_iter = zeros(1,robots);
cases_each_robot_dist = Inf * ones(robots,robots);
tmp_pop_8       = cell(1,8);
new_pop         = cell(1,pop_size);
total_dist      = zeros(1,pop_size);
dist_history    = zeros(1,num_iter);
if show_prog
	show_prog_fig = figure;
    show_prog_fig.Position = [0 0 1000 500] + 50;
end
%==========================================================================
% ----=== TRANSFORMATION --> multiple chromosome [BEGIN] ===----
pop = cell(1,pop_size);
for k = 1: pop_size
	pop{k}.ch{1} = pop_rte(k, 1:pop_brk(k,1));
	for j=2:robots-1
    	pop{k}.ch{j} = pop_rte(k, pop_brk(k,j-1)+1:pop_brk(k,j));
	end
	pop{k}.ch{robots} = pop_rte(k, pop_brk(k,end)+1:n);
end
pop_real_rte = pop; %population of routes including charging points
% ----=== TRANSFORMATION --> multiple chromosome [END] ===----

%% Loop starts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

penalty_rate = 1e2;
start_time = cputime; % get actual time for performance measure
tic;
fprintf('Calculating! \n');
for iter = 1:num_iter  % Evaluate Members of the Population
    total_maxs = zeros(1,pop_size);
    total_mins = inf * ones(1,pop_size);
    pop_real_rte = pop;
    pop_robot_distribution = Inf * ones(pop_size, robots);
    for p = 1:pop_size
        d = 0;
        sman_index = 0;
        for s = 1:length(pop{p}.ch)
            sman = pop_real_rte{p}.ch{s};
            d2 = 0;
            if ~isempty(sman)
                sman_index = sman_index  + 1;
                bateria = max_distance;
                d2 = d2 + cmat(sman(1),1); % Add Start Distance
                bateria = bateria - cmat(sman(1),1);
                if (bateria < 0)
                    d2 = d2 + (d2 - max_distance) * penalty_rate;
                end
                k = 1;
                while (k < length(sman))  %add distance and check in each cycle if it surpassed the distance
                    d2 = d2 + dmat(sman(k),sman(k+1));
                    bateria = bateria - dmat(sman(k),sman(k+1));
                    if (bateria < 0)
                        [d2, bateria, sman] = backtrack_for_battery(k+1,sman);
                        k = k + 1;
                        if (bateria < 0)
                            d2 = d2 + (d2 - max_distance) * penalty_rate;
                        end
                    end
                    k = k + 1;
                end
                % Add End Distance
    %==========================================================================%
    % MULTI DEPOT                                                              %
    %             closest_chrgr = get_closest_charger(sman(end));              %
    %             d2 = d2 + cmat(sman(end),closest_chrgr);                     %
    %             bateria = bateria - cmat(sman(end),closest_chrgr);           %
    %==========================================================================%
    %==========================================================================%
    % SINGLE DEPOT                                                             %
                d2 = d2 + cmat(sman(end),1);                                   %
                bateria = bateria - cmat(sman(end),1);                         %
    %==========================================================================%
                if (bateria < 0 && length(sman) > 1)
                    [d2, bateria, sman] = backtrack_for_battery(length(sman),sman);
                end
    %==========================================================================%
    % MULTI DEPOT                                                              %
    %             sman(end+1)  = -closest_chrgr;                               % 
    %==========================================================================%
    %==========================================================================%
    % SINGLE DEPOT                                                             %
                sman(end+1)  = -1;                                             %
    %==========================================================================%                                                   
                pop_real_rte{p}.ch{s} = sman;                                  
                if (bateria < 0)
                    d2 = d2 + (d2 - max_distance) * penalty_rate;
			    end
			    if d2 > total_maxs(p)
				    total_maxs(p) = d2;
			    end
			    if d2 < total_mins(p)
				    total_mins(p) = d2;
			    end
    
                pop_robot_distribution(p, sman_index) = d2;
            end
            d = d + d2;
        end
        total_dist(p) = d;
    end
    total_dist_F = total_dist;

    %Apply fitness function

    for p = 1:pop_size
	    total_dist_F(p) = (1-weight_wdm)*total_dist_F(p) + weight_wdm*(total_maxs(p) - total_mins(p));
    end

    %% Plot current best solution

    % Find the Best Route in the Population
    [min_dist_F,index] = min(total_dist_F);
    min_dist = total_dist(index);
    dist_history(iter) = min_dist;
    if show_prog % Plot the current best Route
        figure(show_prog_fig);
        subplot(1,2,2)
        opt_rte = pop_real_rte{index}; % the best solution so far
        semilogx(dist_history(1:iter),'b','LineWidth',2);
        set(gca,'XLim',[0 num_iter+1],'YLim',[0 1.1*max([1 dist_history(1:iter)])]);
        if (max(dist_history(1:iter)) >= penalty_rate)
            set(gca, 'YScale', 'log');
        end
        sgtitle(sprintf('MTSP GAC Results History. Iteration = %d', iter))
        title("Total dist / number of iteration ")
        ylabel('Total dist [m]')
        xlabel('Number of iterations [#]')
        subplot(1,2,1)
        plot_results_charge([charging_points(1,:); xy],charging_points(2:end,:),opt_rte.ch,clr)
        title(sprintf('Total Distance = %1.4f, Iteration = %d',min_dist,iter));
    end
    if min_dist_F < global_min
        global_min = min_dist_F; % the optimal solution so far
        opt_rte = pop_real_rte{index}; % the best solution so far
        opt_time = cputime - start_time; % compute the elapsed time
        opt_iter = iter; % store the iteration number
    end
    
    % save best values for each amount of robots
    for i=1:pop_size
	    %first, we get the index by counting the non-empty cells in each member of the population
	    number_of_robot = sum(cellfun(@(x) ~isempty(x), pop{i}.ch)); 
        if total_dist_F(i) < cases_min_F(number_of_robot)
            cases_min(number_of_robot) = total_dist(i); % the optimal solution so far for this case
            cases_min_F(number_of_robot) = total_dist_F(i); % the optimal solution so far for this case
		    cases_routes_min{number_of_robot} = pop_real_rte{i}.ch; % the best solution so far for this case
		    cases_opt_time(number_of_robot) = cputime - start_time; % compute the elapsed time
    	    cases_opt_iter(number_of_robot) = iter; % store the iteration number
            cases_each_robot_dist(number_of_robot, :) = pop_robot_distribution(i,:);
        end 
    end
    
    %% Genetic Algorithm Operators%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rand_grouping = randperm(pop_size);
    for p = 8:8:pop_size
        rpop    = pop(rand_grouping(p-7:p));
        dists   = total_dist(rand_grouping(p-7:p));
        [~,idx] = min(dists);
        best_of_8 = rpop{idx};
        best_of_8.ch(:,cellfun(@(c) isempty(c), best_of_8.ch)) = [];
    
        for k = 1:8 % Generate New Solutions
    
            tmp_pop_8{k} = best_of_8;
            lbestch = length(best_of_8.ch);
            switch k
                case 2 % Flip
                    r = randperm(lbestch);
                    smen = r(1:ceil(rand*lbestch)); % robots selected for flip
                    for k2 = 1:length(smen)
                        if ~isempty(best_of_8.ch{smen(k2)})
                            rte_ins_pts = sort(ceil(length(best_of_8.ch{smen(k2)})*rand(1,2)));
                            I = rte_ins_pts(1);
                            J = rte_ins_pts(2);
                            tmp_pop_8{k}.ch{smen(k2)}(I:J)   = fliplr(tmp_pop_8{k}.ch{smen(k2)}(I:J));
                        end
                    end
                case 3 % Swap
                    smen = ceil(rand(1,2)*lbestch); % the 2 robots selected for swap
                    rte_ins_pts = sort(ceil(min(length(best_of_8.ch{smen(1)}),length(best_of_8.ch{smen(2)}))*rand(1,2)));
                    I = rte_ins_pts(1);
                    J = rte_ins_pts(2);
                    if ~isempty(best_of_8.ch{smen(1)})
                        tempseq = tmp_pop_8{k}.ch{smen(1)}(I:J);
                        tmp_pop_8{k}.ch{smen(1)}(I:J) = tmp_pop_8{k}.ch{smen(2)}(I:J);
                        tmp_pop_8{k}.ch{smen(2)}(I:J) = tempseq;
                    end
                case 4 % Slide
                    r = randperm(lbestch);
                    smen = r(1:ceil(rand*lbestch)); % robots selected for slide
                    toslide = tmp_pop_8{k}.ch{smen(1)}(end);
                    for k2 = 2:length(smen)
                        if ~isempty(best_of_8.ch{smen(k2)})
                            tempgene = tmp_pop_8{k}.ch{smen(k2)}(end);
                            tmp_pop_8{k}.ch{smen(k2)}(2:end) = tmp_pop_8{k}.ch{smen(k2)}(1:end-1);
                            tmp_pop_8{k}.ch{smen(k2)}(1) = toslide;
                            toslide = tempgene;
                        end
                    end
                    tmp_pop_8{k}.ch{smen(1)}(2:end) = tmp_pop_8{k}.ch{smen(1)}(1:end-1);
                    tmp_pop_8{k}.ch{smen(1)}(1) = toslide;
                case 5 % crossover
                    if (lbestch > 1)
                        if (sum(cellfun(@(c) ~isempty(c), tmp_pop_8)) > 1)
                            offsets = crossover_op(tmp_pop_8{k});
                        end
                        tmp_pop_8{k}.ch{offsets{1}(1)} = offsets{2};
                        tmp_pop_8{k}.ch{offsets{1}(2)} = offsets{3};
                    end
                case 6 % Flip and Crossover
                    if (use_complex == 1)
                        r = randperm(lbestch);
                        smen = r(1:ceil(rand*lbestch)); % robots selected for flip
                        for k2 = 1:length(smen)
                            rte_ins_pts = sort(ceil(length(best_of_8.ch{smen(k2)})*rand(1,2)));
                            I = rte_ins_pts(1);
                            J = rte_ins_pts(2);
                            tmp_pop_8{k}.ch{smen(k2)}(I:J)   = fliplr(tmp_pop_8{k}.ch{smen(k2)}(I:J));
                        end
                    end
    
                    % --== CROSSOVER ==--
                    if (lbestch > 1)
                        if (sum(cellfun(@(c) ~isempty(c), tmp_pop_8)) > 1)
                            offsets = crossover_op(tmp_pop_8{k});
                        end
                        tmp_pop_8{k}.ch{offsets{1}(1)} = offsets{2};
                        tmp_pop_8{k}.ch{offsets{1}(2)} = offsets{3};
                    end
                case 7 % Swap and Crossover
                    if (use_complex == 1)
                        smen = ceil(rand(1,2)*lbestch); % the 2 robots selected for swap
                        rte_ins_pts = sort(ceil(min(length(best_of_8.ch{smen(1)}),length(best_of_8.ch{smen(2)}))*rand(1,2)));
                        I = rte_ins_pts(1);
                        J = rte_ins_pts(2);
                        tempseq = tmp_pop_8{k}.ch{smen(1)}(I:J);
                        tmp_pop_8{k}.ch{smen(1)}(I:J) = tmp_pop_8{k}.ch{smen(2)}(I:J);
                        tmp_pop_8{k}.ch{smen(2)}(I:J) = tempseq;
                    end
    
                    % --== CROSSOVER ==--
                    if (lbestch > 1)
                        if (sum(cellfun(@(c) ~isempty(c), tmp_pop_8)) > 1)
                            offsets = crossover_op(tmp_pop_8{k});
                        end
                        tmp_pop_8{k}.ch{offsets{1}(1)} = offsets{2};
                        tmp_pop_8{k}.ch{offsets{1}(2)} = offsets{3};
                    end
                case 8 % Slide and Crossover
                    if (use_complex == 1)
                        r = randperm(lbestch);
                        smen = r(1:ceil(rand*lbestch)); % robots selected for slide
                        toslide = tmp_pop_8{k}.ch{smen(1)}(end);
                        for k2 = 2:length(smen)
                            tempgene = tmp_pop_8{k}.ch{smen(k2)}(end);
                            tmp_pop_8{k}.ch{smen(k2)}(2:end) = tmp_pop_8{k}.ch{smen(k2)}(1:end-1);
                            tmp_pop_8{k}.ch{smen(k2)}(1) = toslide;
                            toslide = tempgene;
                        end
                        tmp_pop_8{k}.ch{smen(1)}(2:end) = tmp_pop_8{k}.ch{smen(1)}(1:end-1);
                        tmp_pop_8{k}.ch{smen(1)}(1) = toslide;
                    end
    
                    % --== CROSSOVER ==--
                    if (lbestch > 1)
                        if (sum(cellfun(@(c) ~isempty(c), tmp_pop_8)) > 1)
                            offsets = crossover_op(tmp_pop_8{k});
                        end
                        tmp_pop_8{k}.ch{offsets{1}(1)} = offsets{2};
                        tmp_pop_8{k}.ch{offsets{1}(2)} = offsets{3};
                    end
            end
        end
        for i=1:8
            new_pop{p-8+i} = tmp_pop_8{i};
        end
    end
    pop = new_pop;
end
toc

%% Final Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if show_res
    for j=1:robots
        figure('Name',"MTSP_GAC | Results with "+j+" robots",'Numbertitle','on');
        if size(charging_points,1) == 1
            plot_results_charge([charging_points; xy],[],cases_routes_min{j},clr)
        else
            plot_results_charge([charging_points(1,:); xy],charging_points(2:end,:),cases_routes_min{j},clr)
        end
        title(sprintf('MTSP GAC Total Distance = %1.1f \n',cases_min(j)));
        daspect([1 1 1])
    end
    figure
    semilogx(dist_history,'b','LineWidth',2);
    title('MTSP GAC Results History: total dist / number of iteration');
    set(gca,'XLim',[0 num_iter+1],'YLim',[0 1.1*max([1 dist_history])]);
    if (max(dist_history(:)) >= 10000)
        set(gca, 'YScale', 'log');
    end
    ylabel('Total dist [m]')
    xlabel('Number of iterations [#]')
end

% Return Outputs
if nargout
    varargout{1} = opt_rte;
    varargout{2} = min_dist;
    varargout{3} = opt_iter;
    varargout{4} = opt_time;
    varargout{5} = dist_history;
	varargout{6} = cases_min;
    varargout{7} = cases_min_F;
	varargout{8} = cases_routes_min;
	varargout{9} = cases_opt_time;
	varargout{10} = cases_opt_iter;
    varargout{11} = cases_each_robot_dist;
end

%% Extra functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %function for executiong the backtrack
    function [distance, battery, modified_route] = backtrack_for_battery(index_panel, route)
        l = index_panel;
        while (l > 1 && route(l) > 0)
            l = l - 1;
        end
        %separates the route into before charging and after charging
        index_chrgr_pnt = l;
        left_rest = route(1:index_chrgr_pnt-1);
        right_rest = route(index_panel+1:end);
        sliced_route = route(index_chrgr_pnt:index_panel);
        if (l == 1)	
            sliced_route = [-1 sliced_route];	
        end
        % compare possibilities of charging
        best_distance = inf;
        best_route = sliced_route;
        for i = 2:length(sliced_route)-1	
            try_route = [sliced_route(1:i) -get_closest_charger(sliced_route(i)) sliced_route(i+1:end)];
            current_distance = get_route_distance(try_route);
            if (current_distance < best_distance)
                best_distance = current_distance;
                best_route = try_route;
            end
        end
        modified_route = [left_rest best_route right_rest];
        if (modified_route(1) < 0)	
            modified_route = modified_route(2:end);	
        end
        [distance, battery] = get_route_distance([-1 modified_route]);
    end

    function [distance, battery] = get_route_distance(route) %get the route distance for every possible path
        distance = 0;
        battery = max_distance;
        for j = 1:length(route)-1
            if (route(j) < 0)
                distance = distance +  cmat(route(j+1),-route(j));
                battery = battery - cmat(route(j+1),-route(j));
            elseif (route(j+1) < 0)
                distance = distance +  cmat(route(j),-route(j+1));
                battery = battery - cmat(route(j),-route(j+1));
                if (battery < 0)
                    distance = distance + (distance - max_distance) * penalty_rate;
                end
                battery = max_distance;
            else
                distance = distance +  dmat(route(j),route(j+1));
                battery = battery - dmat(route(j),route(j+1));
            end
            if (battery < 0)
                distance = distance + (distance - max_distance) * penalty_rate;
            end
        end
    end

    %function for getting the closest charging point to the given panel
    function charger_number = get_closest_charger(panel)
        [~, charger_number] = min(cmat(panel,:));
    end

    % Generate Random Set of Break Points
    function breaks = randbreaks()
        if min_route == 1 % No Constraints on Breaks
            tmp_brks = randperm(n-1);
            breaks = sort(tmp_brks(1:num_brks));
        else % Force Breaks to be at Least the Minimum route Length
            num_adjust = find(rand < cum_prob,1)-1;
            spaces = ceil(num_brks*rand(1,num_adjust));
            adjust = zeros(1,num_brks);
            for kk = 1:num_brks
                adjust(kk) = sum(spaces == kk);
            end
            breaks = min_route*(1:num_brks) + cumsum(adjust);
        end
	end

	% One-point crossover
	function offsets = crossover_op(parent)
		% --== CROSSOVER ==--
		r = randperm(lbestch);
		s_men = r(1:2); % robots selected for crossover
		if sum(cellfun(@length, parent.ch)) < length(xy)-1
			disp('Not enough location!');
		end
		M1 = length(parent.ch{s_men(1)});
		M2 = length(parent.ch{s_men(2)});
		if M1>1
			rr = randperm(M1-1);
		else
			rr = 1;
		end
		cp(1) = rr(1); % point of the crossover in the first robot
		q1 = max(1,cp(1)-M1+min_route); % lower bound of the crossover point for the second robot
		q2 = min(M2,cp(1)+M2-min_route); % upper bound of the crossover point for the second robot
		rr = q1-1+randperm(q2-q1+1);
		cp(2) = rr(1);
		tempseq1 = parent.ch{s_men(1)}(cp(1)+1:end);
		
		offsets = cell(1,3);
		offsets{1} = s_men;
		if (rand <= merging_prob) % Merges the two robots into a single one
			offsets{2} = [parent.ch{s_men(1)} parent.ch{s_men(2)}];
		else
			offsets{2} = [parent.ch{s_men(1)}(1:cp(1)) parent.ch{s_men(2)}(cp(2)+1:end)];
			offsets{3} = [parent.ch{s_men(2)}(1:cp(2)) tempseq1];
		end
    end
end