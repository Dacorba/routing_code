function routes = change_routes_format(routes_cell_array) %function for converting the routes format to be compatible with plots_results
    if isempty(routes_cell_array)
        routes = routes_cell_array;
        return
    end
    robots = length(routes_cell_array);
    routes = {};
    if iscell(routes_cell_array)
        for j = 1:length(routes_cell_array)
            route = routes_cell_array{j};
            if isempty(route)
                continue
            end
            route = route - 1;
            route = [route -1];
            routes{end+1} = route;
        end
    else
        for j = 1:length(routes_cell_array(:,1))
            route = routes_cell_array(j,:);
            if sum(route) == length(route)
                continue
            end
            route = route(route~=1);
            route = route - 1;
            route = [route -1];
            routes{end+1} = route;
        end
    end
end