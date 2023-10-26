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

G = graph(idxs(:,1),idxs(:,2)); %Create a graph where the stops are nodes and the trips are edges
figure
hGraph = plot(G,'XData',PoILon,'YData',PoILat,'LineStyle','none','NodeLabel',{});
hold on

%% Constraints
%Create the linear constraints that each stop has two associated trips

 Aeq = spalloc(PoI,length(idxs),PoI*(PoI-1)); %spalloc Allocate a sparse matrix
for ii = 1:PoI
    whichIdxs = (idxs == ii); %Find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); %Include trips where ii is at either end
    Aeq(ii,:) = whichIdxs'; %Include in the constraint matrix
end
beq = 2*ones(PoI,1);

%% Binary Bounds

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);

%% Optimize Using intlinprog

opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

 x_tsp = logical(round(x_tsp));
Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2),[],numnodes(G));
%Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2)); % Also works in most cases

%% Visualize Solution

hold on
highlight(hGraph,Gsol,'LineStyle','-')
title('Solution with Subtours')

%% Subtour Constraints
%Because you can't add all of the subtour constraints, take an iterative approach

tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % number of subtours
fprintf('# of subtours: %d\n',numtours);
A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];

while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,PoI)]; % A guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1) + 1; % Counter for indexing
        subTourIdx = find(tourIdxs == ii); % Extract the current subtour
        %         The next lines find all of the variables associated with the
        %         particular subtour, then add an inequality constraint to prohibit
        %         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx) - 1; % One less trip than subtour stops
    end
    
    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);
    x_tsp = logical(round(x_tsp));
    Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2),[],numnodes(G));
    % Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2)); % Also works in most cases
    
    % Visualize result
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-')
    drawnow
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % number of subtours
    fprintf('# of robots: %d\n',numtours)
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