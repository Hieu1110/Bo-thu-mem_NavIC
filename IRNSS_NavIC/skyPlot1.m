function hpol = skyPlot1(varargin)

%% Check arguments and sort them ==========================================
[hAxis, args, nargs] = axescheck(varargin{:});

if nargs < 3 || nargs > 4
    error('Requires 3 or 4 data arguments.')
elseif nargs == 3
    [az, el, prn]   = deal(args{1:3});
    line_style      = 'auto';    
else
    [az, el, prn, line_style] = deal(args{1:4});
end

if ischar(az) || ischar(el) || ischar(prn)
    error('AZ and EL must be numeric.');
end

if ~isequal(size(az), size(el))
    error('AZ and EL must be same size.');
end

%% Prepare axis ===========================================================
hAxis = newplot(hAxis);

%--- Get x-axis text color so grid is in same color -----------------------
tc = get(hAxis, 'xcolor');

hold(hAxis, 'on');

%--- Plot white background ------------------------------------------------
rectangle('position', [-90, -90, 180, 180], ...
          'Curvature', [1 1], ...
          'facecolor', 'white', ...
          'edgecolor', tc);

%% Plot spokes ============================================================

%--- Find spoke angles ----------------------------------------------------
% Only 6 lines are needed to divide circle into 12 parts
th = (1:6) * 2*pi / 12;

%--- Convert spoke end point coordinate to Cartesian system ---------------
cst = cos(th); snt = sin(th);
cs = [cst; -cst];
sn = [snt; -snt];

%--- Plot the spoke lines -------------------------------------------------
line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ...
    'handlevisibility', 'off');

%% Annotate spokes in degrees =============================================
rt = 1.1 * 90;

for i = 1:max(size(th))

    %--- Write text in the first half of the plot -------------------------
    text(rt*snt(i), rt*cst(i), int2str(i*30), ...
        'horizontalalignment', 'center', 'handlevisibility', 'off');

    if i == max(size(th))
        loc = int2str(0);
    else
        loc = int2str(180 + i*30);
    end

    %--- Write text in the opposite half of the plot ----------------------
    text(-rt*snt(i), -rt*cst(i), loc, ...
        'handlevisibility', 'off', 'horizontalalignment', 'center');
end

%% Plot elevation grid ====================================================

%--- Define a "unit" radius circle ----------------------------------------
th = 0 : pi/50 : 2*pi;
xunit = cos(th);
yunit = sin(th);

%--- Plot elevation grid lines and tick text ------------------------------
for elevation = 0 : 15 : 90
    elevationSpherical = 90*cos((pi/180) * elevation);

    line(yunit * elevationSpherical, xunit * elevationSpherical, ...
        'lineStyle', ':', 'color', tc, 'linewidth', 0.5, ...
        'handlevisibility', 'off');

    text(0, elevationSpherical, num2str(elevation), ...
        'BackgroundColor', 'white', 'horizontalalignment','center', ...
        'handlevisibility', 'off');
end

%--- Set view to 2-D ------------------------------------------------------
view(0, 90);

%--- Set axis limits ------------------------------------------------------
%save some space for the title
axis([-95 95 -90 101]);

%% Transform elevation angle to a distance to the center of the plot ------
elSpherical = 90*cos(el * pi/180);

%--- Transform data to Cartesian coordinates ------------------------------
yy = elSpherical .* cos(az * pi/180);
xx = elSpherical .* sin(az * pi/180);

%% Define color map =======================================================
% Generate a colormap with as many colors as the number of satellites
num_satellites = size(az, 1);
colors = hsv(num_satellites);  % Use the 'hsv' colormap for more distinct colors

%% Plot data on top of the grid ===========================================

hpol = gobjects(num_satellites, 1);  % Preallocate plot handles
legend_labels = {};  % Initialize empty cell array for legend labels
legend_handles = [];  % Initialize empty array for legend handles

for j = 1:size(az, 2)  % Loop over each time step
    for i = 1:num_satellites
        if xx(i, j) ~= 0 && yy(i, j) ~= 0
            if strcmp(line_style, 'auto')
                %--- Plot with "default" line style -----------------------------------
                hpol(i) = plot(hAxis, xx(i, 1:j)', yy(i, 1:j)', '.-', 'Color', colors(i, :));
            else
                %--- Plot with user specified line style ------------------------------
                hpol(i) = plot(hAxis, xx(i, 1:j)', yy(i, 1:j)', line_style, 'Color', colors(i, :));
            end

            %--- Mark the current position of the satellite -------------------------
            plot(hAxis, xx(i, j)', yy(i, j)', 'o', 'MarkerSize', 7, 'Color', colors(i, :));

            % Add to legend only once per PRN
            if j == 1 && prn(i) ~= 0
                legend_labels{end+1} = ['PRN ', int2str(prn(i))];
                legend_handles(end+1) = hpol(i);
            end
        end
    end
    pause(1);  % Pause for 1 second
end

%--- Place satellite PRN numbers at the latest position -------------------
for i = 1:num_satellites
    if prn(i) ~= 0
        text(xx(i, end), yy(i, end), ['  ', int2str(prn(i))], 'color', colors(i, :));
    end
end

%--- Make sure both axis have the same data aspect ratio ------------------
axis(hAxis, 'equal');

%--- Switch off the standard Cartesian axis -------------------------------
axis(hAxis, 'off');

%--- Add legend to the plot -----------------------------------------------
legend(legend_handles, legend_labels, 'Location', 'bestoutside');

end
