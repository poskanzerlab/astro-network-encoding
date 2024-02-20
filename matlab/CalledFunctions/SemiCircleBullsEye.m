function SemiCircleBullsEye(Struct2Plot, Sig, RANGE, COLORSCHEME,...
    XTICKLABELS, WEDGELABELS, param)
%SemiCircleBullsEye creates a semicircle bullseye plot with each wedge of
%the bullseye colored by a value (magnitude of change frome baseline).
%Area/time bins significantly above baseline will have an asterisk
%   Struct2Plot: a structured array with a field for each band in the
%   bullseye. Each field is a matrix with rows being samples/t-series to
%   average across and columns are time bins (each representing a wedge of
%   the band)
%   Sig: a structured array with the same field names as Struct2Plot. Each
%   field is a 1xn vector of 0s and 1s where n = the same number of columns as matrices
%   in Struct2Plot, determining if that wedge of the band is signficantly above baseline
%   RANGE: an integer to determine the range for colors
%   COLORSCHEME: 2x1 cell array determining the type of color scheme:
%   {'div', 'RdGy'}
%   XTICKLABELS: 
%   WEDGELABELS: 1xn cell array with the strings to label each wedge of the
%   semi-circle
%   param: string describing the parameter being plotted

% =========================================================================================================================================================================
% Michelle Cahill 20220325: calls createBullseye, cbrewer and natsortfiles
% Adrian (2022). Bullseye Plot.zip (https://www.mathworks.com/matlabcentral/fileexchange/47454-bullseye-plot-zip), MATLAB Central File Exchange. Retrieved March 25, 2022.
% =========================================================================================================================================================================

cmap = flipud(cbrewer(COLORSCHEME{1}, COLORSCHEME{2}, RANGE*2,'PCHIP'));
N_AreaBands = length(fieldnames(Struct2Plot));
BAND = natsortfiles(fieldnames(Struct2Plot));

InnerRad = NaN(N_AreaBands, 1);
OuterRad = NaN(N_AreaBands, 1);
nSeg = NaN(N_AreaBands, 1);
startingAngle = NaN(N_AreaBands, 1);
CHANGE = NaN(N_AreaBands, size(Struct2Plot.(BAND{1}),2));
SIG = NaN(N_AreaBands, size(Struct2Plot.(BAND{1}),2));

for band = 1:N_AreaBands
    InnerRad(band) = band - 0.5;
    OuterRad(band) = band + 0.5;
    nSeg(band) = size(Struct2Plot.(BAND{band}), 2);
    startingAngle(band) = -180/nSeg(band);
    CHANGE(band,:) = mean(Struct2Plot.(BAND{band}), 1, 'omitnan'); %first row is near the uncaging site, last row is far
    SIG(band,:) = Sig.(BAND{band});
end
CHANGE_flip = flipud(CHANGE); %first row now corresponds to farthest bin from uncaging and last row is closest bin to uncaging
CHANGE_CAT = reshape(CHANGE_flip', 1,[]);

SIG_flip = flipud(SIG);
SIG_CAT = reshape(SIG_flip', 1,[]);

figure;
set(gcf,'position', [100, 200, 600,500])
BULLSEYE = createBullseye([InnerRad, OuterRad, nSeg, startingAngle]);
set(BULLSEYE,'Color','k','LineWidth',2)
hold on
for b = 1:length(CHANGE_CAT)
    patch(BULLSEYE(b).XData, BULLSEYE(b).YData, cmap(ceil(CHANGE_CAT(b))+RANGE,:));
    if SIG_CAT(b) == 1 %if that bin is signficantly above baseline, plot an asterisk in the middle
        POLYIN = polyshape(BULLSEYE(b).XData, BULLSEYE(b).YData);
        [CENT_X, CENT_Y] = centroid(POLYIN);
        plot(CENT_X, CENT_Y, 'k*')
        clear POLYIN CENT_X CENT_Y
    end
end
clear b
% create a colorbar
colormap(cmap)
c = colorbar('Ticks', [0,0.5,1],...
    'TickLabels', {sprintf('%d', RANGE*-1),'0',sprintf('%d', RANGE)});
c.Label.String = sprintf('change in mean event %s', lower(param));
clear c
% Add in labels for wedges (WEDGELABELS)
for w = 1:length(WEDGELABELS)
    temp_x = BULLSEYE(w).XData(70);
    temp_y = BULLSEYE(w).YData(70);
    if w < length(WEDGELABELS)/2
        addX = -1.25;
        addY = 0.25;
    elseif w == ceil(length(WEDGELABELS)/2)
        addX = -1;
        addY = 0.5;
    else
        addX = 0.1;
        addY = 0.25;
    end
%     text(temp_x, temp_y, 'XX')
    text(temp_x+addX, temp_y+addY, WEDGELABELS{w})
clear temp_x temp_y addX addY
end

% Add in labels for bands in the arc (XTICKLABELS) from outside in

FLIP_XLABELS = flip(XTICKLABELS);
for a = 1:N_AreaBands %length(WEDGELABELS):length(WEDGELABELS):length(CHANGE_CAT)
    band = length(WEDGELABELS) * a;
    inner_x = BULLSEYE(band).XData(1);
    inner_y = BULLSEYE(band).YData(1);
    outer_x = BULLSEYE(band).XData(90);
    outer_y = BULLSEYE(band).YData(90);
    %     if w < length(WEDGELABELS)/2
    %         addX = -1.25;
    %         addY = 0.25;
    %     elseif w == ceil(length(WEDGELABELS)/2)
    %         addX = -1;
    %         addY = 0.5;
    %     else
    %         addX = 0.1;
    %         addY = 0.25;
    %     end
    %     text(inner_x, inner_y, 'X')
    if a == 1
        text(outer_x, outer_y-0.5, num2str(FLIP_XLABELS(a)))
    end
    text(inner_x-0.25, inner_y-0.5, num2str(FLIP_XLABELS(a+1)))
    clear temp_x temp_y addX addY 
end
end