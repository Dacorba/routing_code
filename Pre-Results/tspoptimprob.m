%% Generate random stops

PoI = 100; % define number of PoI
PoILon = round(rand(1,PoI)*100);
PoILat = round(rand(1,PoI)*100); 

%% Calculate Distances Between Points

idxs = nchoosek(1:PoI,2); % Generate all the trips, meaning all pairs of stops
dist = hypot(PoILat(idxs(:,1)) - PoILat(idxs(:,2)), ... % Calculate all the trip distances
             PoILon(idxs(:,1)) - PoILon(idxs(:,2)));
lendist = length(dist);

%% Create Graph and Plot Scenario

G = graph(idxs(:,1),idxs(:,2)); % Create a graph where the stops are nodes and the trips are edges
figure
hGraph = plot(G,'XData',PoILon,'YData',PoILat,'LineStyle','none','NodeLabel',{});
hold on

%% Create Variables and Problem
% Create an optimization problem with binary optimization variables representing the potential trips

tsp = optimproblem;
trips = optimvar('trips',lendist,1,'Type','integer','LowerBound',0,'UpperBound',1);

tsp.Objective = dist*trips;

%% Constraints
% Create the linear constraints that each stop has two associated trips

constr2trips = optimconstr(PoI,1);
for stop = 1:PoI
    whichIdxs = outedges(G,stop); % Identify trips associated with the stop
    constr2trips(stop) = sum(trips(whichIdxs)) == 2;
end

tsp.Constraints.constr2trips = constr2trips;

%% Solve Initial Problem

opts = optimoptions('intlinprog','Display','off');
tspsol = solve(tsp,'options',opts)

%% Visualize Solution

tspsol.trips = logical(round(tspsol.trips));
Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2),[],numnodes(G));
% Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2)); % Also works in most cases

hold on
highlight(hGraph,Gsol,'LineStyle','-')
title('Solution with Subtours')

%% Subtour Constraints
% Because you can't add all of the subtour constraints, take an iterative approach

tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % Number of subtours
fprintf('# of subtours: %d\n',numtours);

k = 1; % Index of added constraints for subtours
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    for ii = 1:numtours
        inSubTour = (tourIdxs == ii); % Edges in current subtour
        a = all(inSubTour(idxs),2); % Complete graph indices with both ends in subtour
        constrname = "subtourconstr" + num2str(k);
        tsp.Constraints.(constrname) = sum(trips(a)) <= (nnz(inSubTour) - 1);
        k = k + 1;        
    end
    
    % Try to optimize again
    [tspsol,fval,exitflag,output] = solve(tsp,'options',opts);
    tspsol.trips = logical(round(tspsol.trips));
    Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2),[],numnodes(G));
    % Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2)); % Also works in most cases
    
    % Plot new solution
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-')
    drawnow

    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % Number of subtours
    fprintf('# of subtours: %d\n',numtours)    
end

pad = 10;
xlim([(min(PoILon)-pad) (max(PoILon)+pad)])
ylim([(min(PoILat)-pad) (max(PoILat)+pad)])
xlabel('X [m]')
ylabel('Y [m]')
hold on
title('TSP');
hold off

%% Solution Quality

disp(output.absolutegap)