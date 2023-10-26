%% Mixed Integer Linear Programing function

function [TotalDist, RobotDist, elapsed_time, robot_trip] = MTSP_MILP(PoI,dist_matrix,robots,clr,PoILon,PoILat,max_distance,weight_wdm,disp_time_bf_milp)

fprintf('Begin MILP \n');
MTSP = optimproblem; 
PoIcoords = [PoILon',PoILat'];
fprintf('Setup done! Computing... \n');
tic;

%% Define Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables are defined using the optivar command. The varaibles are used to 
% describe the problem objective and constraints
% 
% $$x_{\textrm{ijk}} =\left\lbrace \begin{array}{ll}1 & \textrm{Salesman}\;k\;\textrm{passes}\;\textrm{from}\;\textrm{city}\;i\;\textrm{to}\;\textrm{city}\;j\\0 
% & \textrm{else}\end{array}\right.$$

Xijk = optimvar('Xijk', PoI+1,PoI+1,robots, 'Type','integer','LowerBound',0,'UpperBound',1);
z1 = optimvar('z1', 'Type','continuous','LowerBound',0);
z2 = optimvar('z2', 'Type','continuous','LowerBound',0);

%% Define the objective function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use dot notation to define the objective function of the minimisation problem 
% 
% $$f_1 =\sum_{k=1}^m \sum_{i=0}^n \sum_{j=0}^n c_{\textrm{ij}} x_{\textrm{ijk}}$$

F1 = optimexpr(PoI+1,PoI+1,robots);

%% 
% Use for loop to step through each 'page' of the Xijk matrix to complete dot 
% product
% 
% Method > https://uk.mathworks.com/matlabcentral/answers/144400-multiplying-a-2d-matrix-with-each-slice-of-3d-matrix

for k = 1:robots
    F1(:,:,k) = Xijk(:,:,k).*dist_matrix ;
end

%% Fitness Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum all distances in F1 to get total distance travelled by all salesman 
% F = sum(sum(sum(F1)));

F2 = (z1 - z2);
F = (1-weight_wdm)*sum(sum(sum(F1))) + weight_wdm*F2; %fitness calculation

if weight_wdm ~= 0
    z1constraints = optimconstr(robots);
    for k=1:robots
        z1constraints(k) = z1 >= sum(sum(F1(:,:,k)));
    end
    MTSP.Constraints.z1constraints = z1constraints;

    z2constraints = optimconstr(robots);
    for k=1:robots
        z2constraints(k) = z2 <= sum(sum(F1(:,:,k)));
    end
    MTSP.Constraints.z2constraints = z2constraints;
else
end

MTSP.Objective = F;

%% Define Constraints%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Battery constraint, no robot can deplete its battery

batteryConstraints = optimconstr(robots);
for k = 1:robots
    batteryConstraints(k) = sum(sum(Xijk(:,:,k) .* dist_matrix)) <= max_distance;
end
MTSP.Constraints.batteryConstraints = batteryConstraints;

%==========================================================================
% City only visited once
% 
% $$\sum_{k=1}^m \sum_{i=1}^n x_{\textrm{ijk}} =1;\;\;\forall j=1,\ldotp \ldotp 
% \ldotp ,n$$
% 
% $$\sum_{k=1}^m \sum_{j=1}^n x_{\textrm{ijk}} =1;\;\;\forall i=1,\ldotp \ldotp 
% \ldotp ,n$$
% 
XikRows = optimconstr(PoI);
for i = 2:PoI+1
    XikRows(i) = sum(sum(Xijk(i,1:end,:))) == 1;
end

XjkColls = optimconstr(PoI);
for j = 2:PoI+1
    XjkColls(j) = sum(sum(Xijk(1:end,j,:))) == 1;
end

MTSP.Constraints.XikRows = XikRows;
MTSP.Constraints.XjkColls = XjkColls;

% display(XikRows)
% display(XjkColls)

%==========================================================================
%  All salesman start at city 1 
% The sum of row 1 and column 1 for all slices equals the number of salesmen
% 
% $$\sum_{k=1}^m \sum_{i=1}^n x_{\textrm{i0k}} =m$$
% 
% $$\sum_{k=1}^m \sum_{j=1}^n x_{0jk} =m$$

XikStart = sum(sum(Xijk(1,1:end,:))) == robots;
XjkStart = sum(sum(Xijk(1:end,1,:))) == robots;

MTSP.Constraints.XikStart = XikStart;
MTSP.Constraints.XjkStart = XjkStart;

% display(XikStart)
% display(XjkStart)

%==========================================================================
% Each slice (ie salesman) has to have a departing and return leg

XikStartSheet = sum(Xijk(1,:,:)) == 1;
XjkStartSheet = sum(Xijk(:,1,:)) == 1;
% 
MTSP.Constraints.XikStartSheet = XikStartSheet;
MTSP.Constraints.XjkStartSheet = XjkStartSheet;

% =========================================================================
%No subtours
% $\sum_{i\in S} \sum_{j\not\in S} x_{\textrm{ijk}} \ge 1;\;\forall k\in V,\forall 
% S\sqsubseteq C$  - question about if i have done this one!

subpartitions_array = subpartitions(PoI+1);
subpartitionConstraints = optimconstr(length(subpartitions_array));
i = 1;
for subpartition = subpartitions_array
	subpartition = cell2mat(subpartition);
	complement = setdiff(1:PoI+1, subpartition);
	subpartitionConstraints(i) = sum(sum(sum(Xijk(subpartition, complement, :)))) >= 1 ;
	i = i + 1;
end

MTSP.Constraints.subpartitionConstraints = subpartitionConstraints;
% =========================================================================

% show(XikStartSheet)
% show(XjkStartSheet)
% Salesman visiting a city must depart that city
% $$\sum_{i=0}^n x_{\textrm{ipk}} -\sum_{j=0}^n x_{\textrm{pjk}} =0,\;k=1,\ldotp 
% \ldotp \ldotp ,M,\;p=0,\ldotp \ldotp \ldotp ,n$$
% 
ConnectJourney = optimconstr(PoI+1,PoI+1,robots);
for i = 2:PoI+1
        for k = 1:robots
            ConnectJourney(i,j,k) = sum(Xijk(i,:,k)) - sum(Xijk(:,i,k)) == 0;
        end
end
MTSP.Constraints.ConnectJourney = ConnectJourney;

%==========================================================================
% Salesman cannot travel to own city
% Set the diagonal of all matricies to zero

Diagonal  = optimconstr(PoI+1,PoI+1);
for i = 1:PoI+1
    for k = 1:robots 
    Diagonal(i,i,k) = Xijk(i,i,k) == 0;
    end
end
MTSP.Constraints.Diagonal = Diagonal;

%% Solve the Problem
% Use the solve command to solve the optimisation platform 
% try
sol = solve(MTSP);

%% Display the results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show the results of the optimisation problem
% 
% Create matrix of the distances travelled

Dist = dist_matrix.*sol.Xijk;

%% 
% Calculate total distance travelled by all salesman

TotalDist = sum(sum(sum(Dist)))

%% 
% Calculate distance travelled by each salesman

RobotDist = zeros(robots,1) ;
for k = 1:robots
    RobotDist(k) = sum(sum(sum(sum((dist_matrix.*sol.Xijk(:,:,k))))));
end

RobotDist

elapsed_time = toc;

%% Plot results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if disp_time_bf_milp == 0
    PoI = PoI+1;
    plot_results_no_charge(PoIcoords,[])
    hold on
    robot_trip=ones(robots,1);
    solXijk = discretize(sol.Xijk(:,:,:), [0.9,1.1]) == 1; % copy the solution, but only showing values close to 1
    for k=1:robots
        i = 1;
        y = 1;
        while sum(sum(solXijk(:,:,k))) ~= 0
            j = 2;
            while j <= PoI && ~(i == 1 && j == 1)
                if (solXijk(i, j, k) ~= 0)
                    y = y + 1;
                    solXijk(i, j, k) = 0;
                    robot_trip(k,y) = j;
                    hold on
                    plot( [PoILon(i) PoILon(j)], [PoILat(i) PoILat(j)],'LineWidth',1.5 , 'Color', clr(k,:))
                    i = j;
                    j = 1;
                    continue;
                end
                j = j + 1;
            end
            i = i + 1;
        end
    end
    title(sprintf('Total Distance = %1.1f',TotalDist));
    daspect([1 1 1])
end

%% Extra function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function result_array = subpartitions(N) %function that creates all combinations until N
        result_array = {};
        for i = 1:(N-1)
            i_cases = nchoosek(1:N, i);
            for j = 1:length(i_cases)
                result_array{end+1} = i_cases(j,:);
            end
        end
    end

end