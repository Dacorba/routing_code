function routes = change_routes_format2(routes)% function to add the depot at the start and the end of the route
    if isempty(routes)
        routes = routes;
        return
    end
    for i = 1:length(routes)
        routes{i} = [0 routes{i}];
        routes{i}(end) = 0;
    end
end

