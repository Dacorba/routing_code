function [] = plot_points_scenario(PoILon,PoILat) %function to plot the scenario only

scatter(PoILon(1),PoILat(1),200,'black',"filled"); %plot depot
hold on
scatter(PoILon(2:end),PoILat(2:end),200,'blue',"square","filled"); %plot PoI
hold on

%% Numeration

text(PoILon(1)+5, PoILat(1), num2str(0), ... %adding depot identification
'Color', 'black', 'FontSize',14)
hold on

PoI = length(PoILon);
for i=2:PoI %adding panel identification
    if i < 11
        text(PoILon(i)-7, PoILat(i), num2str(i-1), ...
        'Color', 'blue', 'FontSize',14)
        hold on
    elseif i > 10 && i < 101
        text(PoILon(i)-10, PoILat(i), num2str(i-1), ...
        'Color', 'blue', 'FontSize',14)
        hold on
    else
        text(PoILon(i)-14, PoILat(i), num2str(i-1), ...
        'Color', 'blue', 'FontSize',14)
        hold on
    end
end
    
%% Graphics Axis Limits
pad = 10;
xlim([(min(PoILon)-pad) (max(PoILon)+pad)])
ylim([(min(PoILat)-pad) (max(PoILat)+pad)])
xlabel('X [m]')
ylabel('Y [m]')
hold on
daspect([1 1 1])

end		    