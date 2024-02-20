function ret = aqua_extract( inputPath, options )
% AQUA_EXTRACT  Extract AQuA events from a res file; called by AQUA_EXPORT
%   ret = AQUA_EXTRACT(inputPath) loads the res-file at `inputPath` and
%           puts all of the data for export into the structure `ret`
%
%   TODO This doesn't handle the no-uncaging case without manual edits; add.
%
%   See also AQUA_EXPORT.


    arguments
        inputPath {mustBeFile}

        options.metadata (1,1) struct = struct
    end


    %% Load data
    %  ---------
    res = load( inputPath, 'res' ).res;


    %% Determine global parameters
    %  ---------------------------

    nCells = size( res.ftsFilter.region.cell.incluLmk, 1 );
    nEvents = size( res.ftsFilter.loc.t0, 2 );

    frameT = res.opts.frameRate;
    pixelUM = res.opts.spatialRes;
    imageSize = res.opts.sz; % in *pixels*

    colCellNames = res.ftsFilter.region.cell.name; % these are strings

    uncageCellCol = find( res.ftsFilter.region.cell.incluLmk == 1 );
    uncageCell = str2num( colCellNames{uncageCellCol} );
    uncagePixelX = res.ftsFilter.region.landMark.center(1);
    uncagePixelY = res.ftsFilter.region.landMark.center(2);

    % Decode direction order
    directionOrder = res.ftsFilter.notes.propDirectionOrder;
    for iDir = 1:4
        if strcmp( directionOrder{iDir}, 'Anterior' )
            anteriorIndex = iDir;
        end
        if strcmp( directionOrder{iDir}, 'Posterior' )
            posteriorIndex = iDir;
        end
        if strcmp( directionOrder{iDir}, 'Left' )
            leftIndex = iDir;
        end
        if strcmp( directionOrder{iDir}, 'Right' )
            rightIndex = iDir;
        end
    end


    %% Determine event locations and membership
    %  ----------------------------------------

    % Onset time for all events
    frameEventAll = res.ftsFilter.loc.t0;

    % Cell identity for each event

    % % Original method
    % cellEventAll = zeros( nEvents, 1 );
    % for i = 1:nEvents
    %     cellEventAll(i) = find( res.ftsFilter.region.cell.memberIdx(i, :) == 1 );
    % end

    % Updated method with sub-selections
    % TODO
    cellEventAll = zeros( nEvents, 1 );
    for i = 1:nEvents
        cellColCur = find( res.ftsFilter.region.cell.memberIdx(i, :) == 1 );
        cellEventAll(i) = str2num( colCellNames{cellColCur} );
    end


    %% Collate event marks
    %  -------------------

    % Marks for each event
    marks = struct;

    marks.area = res.ftsFilter.basic.area';
    marks.peri = res.ftsFilter.basic.peri';
    marks.circMetric = res.ftsFilter.basic.circMetric';

    marks.dffMax = res.ftsFilter.curve.dffMax';
    marks.dffMax2 = res.ftsFilter.curve.dffMax2';
    marks.rise19 = res.ftsFilter.curve.rise19';
    marks.fall91 = res.ftsFilter.curve.fall91';
    marks.width55 = res.ftsFilter.curve.width55';
    marks.width11 = res.ftsFilter.curve.width11';
    marks.decayTau = res.ftsFilter.curve.decayTau';

    marks.propGrowAnterior = res.ftsFilter.propagation.propGrowOverall(:, anteriorIndex);
    marks.propGrowPosterior = res.ftsFilter.propagation.propGrowOverall(:, posteriorIndex);
    marks.propGrowLeft = res.ftsFilter.propagation.propGrowOverall(:, leftIndex);
    marks.propGrowRight = res.ftsFilter.propagation.propGrowOverall(:, rightIndex);
    marks.propShrinkAnterior = res.ftsFilter.propagation.propShrinkOverall(:, anteriorIndex);
    marks.propShrinkPosterior = res.ftsFilter.propagation.propShrinkOverall(:, posteriorIndex);
    marks.propShrinkLeft = res.ftsFilter.propagation.propShrinkOverall(:, leftIndex);
    marks.propShrinkRight = res.ftsFilter.propagation.propShrinkOverall(:, rightIndex);

    marks.nOccurSameTime = res.ftsFilter.networkAll.nOccurSameTime;
    marks.nOccurSameLoc = res.ftsFilter.networkAll.nOccurSameLoc(:, 1);
    marks.nOccurSameLocSize = res.ftsFilter.networkAll.nOccurSameLoc(:, 2);

    % TODO Multiple landmarks?
    marks.uncageDistMin = res.ftsFilter.region.landmarkDist.distMin;
    marks.uncageDistAvg = res.ftsFilter.region.landmarkDist.distAvg;

    % Centroid
    marks.centroidXPixels = zeros( size( marks.area ) );
    marks.centroidYPixels = zeros( size( marks.area ) );
    marks.centroidXUM = zeros( size( marks.area ) );
    marks.centroidYUM = zeros( size( marks.area ) );
    marks.centroidRelXUM = zeros( size( marks.uncageDistMin ) );
    marks.centroidRelYUM = zeros( size( marks.uncageDistMin ) );
    for i = 1 : numel( marks.centroidXPixels )
        curX2D = res.ftsFilter.loc.x2D{1, i};
        [r, c] = ind2sub( imageSize(1:2), curX2D );
        startX = min( c );
        startY = min( r );
        
        curMap = res.ftsFilter.basic.map{1, i};
        relCentroid = regionprops( true( size( curMap ) ), curMap, 'WeightedCentroid').WeightedCentroid;
        centroid = [startX, startY] + relCentroid;
        
        marks.centroidXPixels(i) = centroid(1);
        marks.centroidYPixels(i) = centroid(2);
        marks.centroidXUM(i) = centroid(1) * pixelUM;
        marks.centroidYUM(i) = centroid(2) * pixelUM;
        marks.centroidRelXUM(i) = (centroid(1) - uncagePixelX) * pixelUM;
        marks.centroidRelYUM(i) = (centroid(2) - uncagePixelY) * pixelUM;
    end
        
    %% Collect exported data
    %  ---------------------
    ret = struct;
        
    ret.fs = 1 / frameT;
    ret.pixelUM = pixelUM;
    ret.imageSize = imageSize;

    ret.uncageCell = uncageCell;
    % Old method, now covered with metadata argument
    % ret.uncageFrameStart = options.uncageFrameStart;
    % ret.uncageFrameEnd = options.uncageFrameEnd;
    ret.uncagePixelX = uncagePixelX;
    ret.uncagePixelY = uncagePixelY;

    ret.eventFrames = frameEventAll';
    ret.eventCells = cellEventAll;

    % Save metadata in a systematic way
    metadataNames = fieldnames( options.metadata );
    for kMeta = 1 : numel( metadataNames )
        metadataData = options.metadata.(metadataNames{kMeta});
        metadataKey = metadataNames{kMeta};
        ret.(metadataKey) = metadataData;
    end

    % Save marks in systematic way
    markNames = fieldnames( marks );
    for kMark = 1 : numel( markNames )
        markData = marks.(markNames{kMark});
        markKey = ['mark_', markNames{kMark}];
        ret.(markKey) = markData;
    end

end
%