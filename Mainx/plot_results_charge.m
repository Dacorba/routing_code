function [] = plot_results_charge(xy,extra_charging_points,routes,clr) %function to plot the results of charging points algorithms results

extra_charging_points = [xy(1,1:2); extra_charging_points];
xy = xy(2:end,1:2);

panels = size(xy, 1); %ALFREDO
PoILon = xy(:,1);
PoILat = xy(:,2);

%plotting PoI
scatter(PoILon, PoILat,200,'blue',"square","filled")
xlabel('X [m]')
ylabel('Y [m]')
hold on

scatter(extra_charging_points(:,1),extra_charging_points(:,2),200,'red',"filled") %plot charging points
%adding charging point identification
text(extra_charging_points(1,1), extra_charging_points(1,2)-6, num2str(0),'Color', 'red', 'FontSize', 12)
for z = 2:size(extra_charging_points,1)
    text(extra_charging_points(z,1), extra_charging_points(z,2)-6, num2str(-(z-1)), ...
        'Color', 'red', 'FontSize', 12)
end


% the graphics axis limits
pad = 10;
xvalues = [PoILon; extra_charging_points(:,1)];
xvalues = xvalues(~isinf(xvalues));
yvalues = [PoILat; extra_charging_points(:,2)];
yvalues = yvalues(~isinf(yvalues));

xlim([(min(xvalues)-pad) (max(xvalues)+pad)])
ylim([(min(yvalues)-pad) (max(yvalues)+pad)])
xlabel('X [m]')
ylabel('Y [m]')
hold on
  
for j = 1:length(routes)
    path = routes{j};
    if isempty(path)
        continue
    end
    if(path(end) ~= -1)
        path = path - 1;
        path = [path -1];
    end
    complete_path = [-1 path];
    for i = 1:length(complete_path)-1
        cdx = zeros(1,2);
        cdy = zeros(1,2);
        for z = [0 1]
            if (complete_path(i+z) > 0)
               cdx(z+1) = PoILon(complete_path(i+z));
               cdy(z+1) = PoILat(complete_path(i+z));
            else
                cdx(z+1) = extra_charging_points(-complete_path(i+z),1);
                cdy(z+1) = extra_charging_points(-complete_path(i+z),2);
            end
        end
        plot(cdx, cdy, 'LineWidth',1.5 , 'Color', clr(j,:))
    end    
end   
end