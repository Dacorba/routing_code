%% MAIN CODE
% The main code allows to: clean work space, specify variables, start the scenario,
% calculate a solution for Brute Force, Random, Mixed Integer Linear Programming and GAs
%
% Important notes to consider:
%
% Note1: BF, Random, MILP and 1 GA do not use charging points, while 2 GA's use them
%
% Note2: The depot location is taken to be the first PoI point
%
% Note3: Deterministic algorithms will find the best solution but demanding 
% computational times (BF and MILP) while the heuristics Random and GA will try to find 
% a near best solution in a computational time controled by the number of iterations, 
% decided by the user
%
% Note4: According to the algorithm, functions called "change_routes" adjust the 
% routes so all algorithms can represent equally their routing solutions
%
% Note5: In the charging points algorithms, there is an option to test multi
% depot, turning every charging points into a possible depot (altought, keeping the first 
% depot as the starting one). To achieve this, the user needs to comment and uncomment 
% specific code lines in those algorithms. These codes lines are sinalized in both 
% scripts: MTSP_GAC_MSD and MTSP_GADC_MSD; the MSD means Multiple/Single Depot
%
%
% Information about the algorithms
%
% Scenary function - the scenario is chosen according to the PoI number;
% there are 6 scenarios previously done to be tested
%
% Plot functions - there are three functions available to plot the 
% different algorithm's solutions : only to plot the scenario, the results
% without charging points, and the results with charging points.
%
% MTSP BF - BF runs once only, because calculates every robot combination at once
%
% MTSP Random - Random runs the same times as the robots' number, plus
% will do this for the number of times choosen by user in the inputs setup
%
% MTSP MILP - MILP runs the same times as the robots' number
%
% MTSP GA - GA runs the same times as the robots' number -1, because
% it needs at least 2 robots to run, opposing to the previous algorithms; 
% although it can find solutions for 1 robot
%
% MTSP GAC - similar to the GA, but includes fixed charging points, whose
% locations are to be defined by the user
%
% MTSP GADC - similar to the GA, but includes charging points, whose
% locations are to be defined by the algorithm, which will select the
% charging points locations dynamically 
%
%
% Author: Daniel Correia Batista
% Email: daniel.batista@tecnico.ulisboa.pt
% Release Date: 28/10/2023
% *************************************************************************
% --== Reference notice ==--
% If you use this implementation in your work, please cite this paper:
%
% Daniel Correia Batista: Routing problem optimization for mobile robots 
% including the number and location of charging points, 2023.
% *************************************************************************
%% Cleaning workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Inputs setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overal inputs============================================================
PoI = 101; % number of points of interest (PoI-1), with depot included
robots = 4; % total number of robots available
MATRIX = cell(robots,6); % MATRIX that saves the outputs for the 6 algorithms
random = 0; % 0 if false, 1 if true, for random scenary
max_distance = 1500; % battery limit
weight_wdm = 0; % determines how much importance is given in the balance of 
%               minimizing the distance between the greatest and smallest
%               routes, and minimizing the total distance
format shortG
%==========================================================================
% Iteration/time based calculation=========================================
time =  0; % 0 based on iteration, 1 based on time
if time == 1
    max_time = 10; % define max time to run algorithms
else
    max_time = 0;
end
%==========================================================================
% Display the computational time taken from 3 to PoI: applied in BF and MILP
disp_time_bf_milp = 1; % display if 1,dont display if 0
%==========================================================================
% Random inputs============================================================
number_of_random_iterations = 2; % number of times random will run for r amount of robots
indiv_results = 0; % 1 to show individual results, 0 to not
indiv_hist = 0; % 1 to show individual evolution for the current amount of robots, 0 to not
if PoI == 6 % define iterations
    iter_tests = [1 10 100 1000 10000];
elseif PoI == 21
    iter_tests = [1 10 100 1000 10000 100000 1000000 10000000]; 
elseif PoI == 51
    iter_tests = [1 10 100 1000 10000 100000 1000000 10000000];
elseif PoI == 101
    iter_tests = [1 10 100 1000 10000 100000 1000000 10000000];
else
    iter_tests = [1 10 100 1000 10000 100000];
end
if time == 1 %if time, define inf iterations to run on time
	iter_tests = inf;
end
%==========================================================================
% GA inputs================================================================
pop_size = 80; % number of solutions per population, must be multiple of 8!!!
num_iter = 200; % number of iterations
num_iter_time = inf; % 
use_complex = 0; % is the flag wether to use complex mutation operators or not
show_prog = 1; % shows the GA progress if true
show_res = 0; % shows the GA results if true
%==========================================================================
% Charging points inputs===================================================
additional_charging_points = []; %vector with coordinates for MTSP GAC
amount_of_new_chargers = 2; % number of available chargers for MTSP GADC
%==========================================================================
% Scenario + dist_matrix===================================================
[PoILon, PoILat] = scenary(PoI, random); % scenario function
PoIcoords = [PoILon;PoILat];
dist_matrix = zeros(PoI-1,PoI-1);
for p = 1:PoI
    for l=p+1:PoI
        dist_matrix(p,l) = norm(PoIcoords(:,p)-PoIcoords(:,l));
        dist_matrix(l,p) = dist_matrix(p,l);
    end
end
%==========================================================================
% Colour===================================================================
clr = [ 1 0.7 0; 0 1 0; 0 1 1; 1 0 1;];
while size(clr,1) < robots
    clr(end+1,:) = [rand rand rand];
end
%==========================================================================

%% Plot Scenario %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot_points_scenario(PoILon, PoILat)
daspect([1 1 1])

%% BF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PoI < 10
% Results for PoI
    [min_dist, best_trip_per_robot_cases, time_elapsed, min_dist_cases] = MTSP_BF(PoI,dist_matrix,robots,clr,PoILon,PoILat,max_distance,weight_wdm,0);
    brute_force_total_time = toc;
    min_distbf = min_dist;
    best_trip_per_robot_casesbf = best_trip_per_robot_cases;
    time_elapsedbf = time_elapsed;
    for i = 1: robots
        MATRIX{i,1} = {min_dist_cases(i) [zeros(robots, 1) best_trip_per_robot_cases(:,:,i)-1] time_elapsed(i)};
    end
% Time taken to calculate from 3 to PoI 
    if disp_time_bf_milp == 1
        figure();
        timesbf = [];
        pcountbf = [];
        for i = 3:PoI
            [min_dist, best_trip_per_robot_cases, time_elapsed, min_dist_cases] = MTSP_BF(i,dist_matrix,robots,clr,PoILon,PoILat,inf,0,disp_time_bf_milp);
            timesbf(end+1)=sum(time_elapsed);
            pcountbf(end+1)=i;
        end
        timesbf
        bar(pcountbf, timesbf);
        set(gca, 'YScale', 'log');
        xlim([2.5 (PoI+0.5)]);
        xlabel('Number of Panels [#]');
        ylabel('Time [T]');
        sgtitle('BF - Time vs. Number of Panels')
    end
end

%% Random %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_best_distance = Inf * ones(1,robots);
current_best_distance_robots = Inf * ones(robots);
current_best_trip_per_robot = ones(robots, PoI + 1, robots);

for j = 1:number_of_random_iterations
    for i = 1:robots
        try
            [min_dist_rand, best_trip_per_robot_rand, best_dist_travelled_robot, time_elapsed_rand, iter_tests] = MTSP_RANDOM(PoI,dist_matrix,i, robots,clr,PoILon,PoILat,max_distance,weight_wdm,time,max_time,iter_tests,indiv_results,indiv_hist);
            if min_dist_rand < current_best_distance(i)
                current_best_trip_per_robot(1:i,:,i) = best_trip_per_robot_rand;
                current_best_distance(i) = min_dist_rand;
                current_best_distance_robots(i,1:i) = best_dist_travelled_robot;
                MATRIX{i,2} = {min_dist_rand, change_routes_format2(change_routes_format(best_trip_per_robot_rand)), time_elapsed_rand};
            end
        catch E
            if ~(E.message == "Unrecognized function or variable 'robot_routes'.")
                errorMessage = sprintf('Error in function %s() at line %d.\n', ...
            E.stack(1).name, E.stack(1).line);
                fprintf(2, '%s\n', errorMessage);
                throw(E)
            end
        end
    end
    figure(105)
    sgtitle('MTSP Random results: total dist / number of iteration')
end

f=figure;
sgtitle('MTSP Random');

for i = 1:robots
    subplot(2,ceil(robots/2),i);
    plot_results_no_charge(PoIcoords',change_routes_format(current_best_trip_per_robot(1:i,:,i)),clr);
    title(sprintf('Total Distance = %1.1f\n',current_best_distance(i)))
    daspect([1 1 1])
end

sprintf("Best distances for each robot:")
display(current_best_distance_robots)

%% MILP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PoI < 17
% Results for PoI
f=figure;
    for i = 1:robots
        subplot(2,ceil(robots/2),i);
        try
            [TotalDist_solver, robotdist, elapsed_time_solver, best_trip_per_robot_solver] = MTSP_MILP(PoI-1,dist_matrix,i,clr,PoILon,PoILat,max_distance,weight_wdm,0);
            best_trip_per_robot_solver = best_trip_per_robot_solver-1;
            best_trip_per_robot_solver(best_trip_per_robot_solver == -1) = 0;
            MATRIX{i,3} = {TotalDist_solver, best_trip_per_robot_solver, elapsed_time_solver};
        catch E
                if (E.message == "Arrays have incompatible sizes for this operation.")              
                    % plot_points_scenario(PoILon, PoILat)   
                    plot_results_no_charge(PoIcoords',{},clr)
                    title(sprintf('Total Distance = %1.1f',Inf));
                    daspect([1 1 1])
                else
                    errorMessage = sprintf('Error in function %s() at line %d.\n', ...
                E.stack(1).name, E.stack(1).line);
                    fprintf(2, '%s\n', errorMessage);
                    throw(E)
                end
        end
    end
    sgtitle('MTSP MILP - Mixed Integer Linear Programming')
%Time taken to calculate from 3 to PoI 
    if disp_time_bf_milp == 1
        figure();
        timesmilp= [];
        pcount = [];
        for i = 3:PoI
            try
                [TotalDist_solver, robotdist, elapsed_time_solver] = MTSP_MILP(i-1,dist_matrix(1:i,1:i),1,clr,PoILon,PoILat,1e20,0,disp_time_bf_milp);
                timesmilp(end+1)=elapsed_time_solver;
            catch E
                if (E.message == "Arrays have incompatible sizes for this operation.")
                    timesmilp(end+1)=0;
                else
                    display("unhandled error")
                    errorMessage = sprintf('Error in function %s() at line %d.\n', ...
                    E.stack(1).name, E.stack(1).line);
                    fprintf(2, '%s\n', errorMessage);
                    throw(E)
                end
            end
            pcount(end+1)=i;
        end

        bar(pcount, timesmilp);
        set(gca, 'YScale', 'log');        
        xlim([2.5 (PoI+0.5)]);
        xlabel('Number of Panels [#]');
        ylabel('Time [T]');        
        sgtitle('MILP - Time vs. Number of Panels')

    end
end

%% MTSP GA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_best_min_dist = inf * ones(1,robots);
current_best_min_dist_F = inf * ones(1,robots);
current_best_opt_route = cell(1,robots);
current_best_cpu_time = zeros(robots,1);
current_best_robot_dist = Inf * ones(robots);

for i=2:robots
    if time == 0
        [opt_rte, min_dist, opt_iter, ...
            opt_time, dist_history, cases_min_dist, cases_min_dist_F, cases_routes_min, ...
            cases_opt_time, cases_opt_ite, cases_each_robot_dist] ...
        = MTSP_GA_SD(PoIcoords',dist_matrix,i,1,max_distance,weight_wdm,pop_size,num_iter,use_complex,show_prog,show_res,clr);      
        for j = 1: length(cases_min_dist)
            if cases_min_dist_F(j) < current_best_min_dist_F(j)  
                current_best_min_dist_F(j) = cases_min_dist_F(j);
                current_best_min_dist(j) = cases_min_dist(j);
                current_best_opt_route{j} = cases_routes_min{j};
                current_best_cpu_time(j) = cases_opt_time(j);
                current_best_robot_dist(j,1:j) = cases_each_robot_dist(j,1:j);
            end
        end
    else
        [opt_rte, min_dist, opt_iter, ...
            opt_time, dist_history, cases_min_dist, cases_min_dist_F, cases_routes_min, ...
            cases_opt_time, cases_opt_ite, cases_each_robot_dist] ...
        = MTSP_GA_SD_time(PoIcoords',dist_matrix,i,1,max_distance,weight_wdm,max_time,pop_size,num_iter_time,use_complex,show_prog,show_res,clr);
        for j = 1: length(cases_min_dist)
            if cases_min_dist_F(j) < current_best_min_dist_F(j) &&  cases_min_dist_F(j) < 1e9
                current_best_min_dist_F(j) = cases_min_dist_F(j);
                current_best_min_dist(j) = cases_min_dist(j);
                current_best_opt_route{j} = cases_routes_min{j};
                current_best_robot_dist(j,1:j) = cases_each_robot_dist(j,1:j);
            end
        end
    end
end

for i = 1:robots
    MATRIX{i,4} = {current_best_min_dist(i), change_routes_format2(change_routes_format(current_best_opt_route{i})),current_best_cpu_time(i)};
end

f=figure;
for j = 1:robots
    subplot(2,ceil(robots/2),j);
    plot_results_no_charge(PoIcoords',current_best_opt_route{j},clr)
    title(sprintf('Total Distance = %1.1f \n',current_best_min_dist(j)));
    daspect([1 1 1])
end

sprintf("Best distances for each robot:")
current_best_robot_dist
sgtitle('MTSP GA')

 %% MTSP GAC + Multi/Single Depot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_best_min_dist = inf * ones(1,robots);
current_best_min_dist_F = inf * ones(1,robots);
current_best_opt_route = cell(1,robots);
current_best_cpu_time = zeros(robots,1);
current_best_robot_dist = Inf * ones(robots);

for i=2:robots
    [opt_rte, min_dist, opt_iter, ...   
    opt_time, dist_history, cases_min_dist, cases_min_dist_F, cases_routes_min, ...
    cases_opt_time, cases_opt_ite, cases_each_robot_dist] ...
    = MTSP_GAC_MSD(PoIcoords',dist_matrix(2:end,2:end),additional_charging_points,i,1,max_distance,weight_wdm,pop_size,num_iter,use_complex,show_prog,show_res,clr);
    for j = 1:length(cases_min_dist)
        if cases_min_dist_F(j) < current_best_min_dist_F(j)  
            current_best_min_dist_F(j) = cases_min_dist_F(j);
            current_best_min_dist(j) = cases_min_dist(j);
            current_best_opt_route{j} = cases_routes_min{j};
            current_best_cpu_time(j) = cases_opt_time(j);
            current_best_robot_dist(j,1:j) = cases_each_robot_dist(j,1:j);
        end
    end
end

f=figure;
for j = 1:robots
    subplot(2,ceil(robots/2),j);
    plot_results_charge(PoIcoords',additional_charging_points,current_best_opt_route{j},clr)
    title(sprintf('Total Distance = %1.1f \n',current_best_min_dist(j)));
    daspect([1 1 1])
end

for i = 1:robots
    MATRIX{i,5} = {current_best_min_dist(i), change_routes_format2(current_best_opt_route{i}),current_best_cpu_time(i)};
end

sprintf("Best distances for each robot:")
current_best_robot_dist
sgtitle('MTSP GAC')


%% MTSP GADC + Multi/Single Depot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_best_min_dist = inf * ones(1,robots);
current_best_min_dist_F = inf * ones(1,robots);
current_best_opt_route = cell(1,robots);
current_best_new_chargers = inf*ones(robots,2,amount_of_new_chargers);
current_best_cpu_time = zeros(robots,1);
current_best_robot_dist = Inf * ones(robots);

for i=4:robots
    [opt_rte, min_dist, opt_iter, ...   
    opt_time, dist_history, cases_min_dist, cases_min_dist_F, cases_routes_min, ...
    cases_opt_time, cases_opt_ite, cases_new_charger, cases_each_robot_dist] ...
    = MTSP_GADC_MSD(PoIcoords',dist_matrix(2:end,2:end),amount_of_new_chargers,i,1,max_distance,weight_wdm,pop_size,num_iter,use_complex,show_prog,show_res,clr);
    for j = 1:length(cases_min_dist)
        if cases_min_dist_F(j) < current_best_min_dist_F(j)  
            current_best_min_dist_F(j) = cases_min_dist_F(j);
            current_best_min_dist(j) = cases_min_dist(j);
            current_best_opt_route{j} = cases_routes_min{j};
            current_best_new_chargers(j,:,:) = cases_new_charger(j,:,:);
            current_best_cpu_time(j) = cases_opt_time(j);
            current_best_robot_dist(j,1:j) = cases_each_robot_dist(j,1:j);
        end
    end
end

f=figure;
for j = 1:robots
    subplot(2,ceil(robots/2),j);
    plot_results_charge(PoIcoords',squeeze(current_best_new_chargers(j,:,:))',current_best_opt_route{j},clr)
    title(sprintf('Total Distance = %1.1f \n',current_best_min_dist(j)));
    daspect([1 1 1])
end

for i = 1:robots
    MATRIX{i,6} = {current_best_min_dist(i), change_routes_format2(current_best_opt_route{i}),current_best_cpu_time(i)};
end

sprintf("Best distances for each robot:")
current_best_robot_dist
sgtitle('MTSP GADC')