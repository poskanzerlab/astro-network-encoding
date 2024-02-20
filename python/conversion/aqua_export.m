function ret = aqua_export( options )
% AQUA_EXPORT  Extract and save AQuA events from a dataset
%   ret = AQUA_EXPORT(options)
%
%   Available keyword options:
%   * 'inputPath': Path to the folder with input aqua data
%   * 'outputPath': Path to the folder to place results
%   * 'experiment': Presets for uncaging start/end frames used in the paper:
%       * 'uncaging': start 215, end 217
%       * 'repeat-25': start 164, end 166
%       * 'repeat-70': start 213, end 215
%   * 'uncageFrameStart', 'uncageFrameEnd': If 'experiment' isn't set, allows
%       manual specification of uncaging start and end frames
%   * 'saveSuffix': Suffix to append to exported event files
%   * 'workers': Number of parallel workers to use (note: each worker can use
%       quite a bit of memory)
%
%   TODO This doesn't handle the no-uncaging case without manual edits; add.
%
%   See also AQUA_EXTRACT.

    arguments
        options.inputPath {mustBeFolder}
        options.outputPath (1,1) string = fullfile( '.', 'export' )
        options.experiment (1,1) string = 'uncaging'
        options.uncageFrameStart (1,1) int = 215
        options.uncageFrameEnd (1,1) int = 217
        options.saveSuffix (1,1) string = 'events'
        options.workers (1,1) int8 = 1
    end


    % Generate metadata
    % TODO This hard-codes the start and end frames; 
    metadata = struct;
    if strcmp( options.experiment, 'uncaging' )
        metadata.uncageFrameStart = 215;
        metadata.uncageFrameEnd = 217;
    elseif strcmp( options.experiment, 'repeat-25' )
        metadata.uncageFrameStart = 164;
        metadata.uncageFrameEnd = 166;
    elseif strcmp( options.experiment, 'repeat-70' )
        metadata.uncageFrameStart = 213;
        metadata.uncageFrameEnd = 215;
    else
        metadata.uncageFrameStart = options.uncageFrameStart;
        metadata.uncageFrameEnd = options.uncageFrameEnd;
    end

    % Make sure the save directory exists
    mkdir( options.outputPath );

    % Enumerate subdirs of the inputPath
    inputPathContents = dir( options.inputPath );
    inputPathDirFlags = [inputPathContents.isdir];
    inputPathSubdirs = inputPathContents( inputPathDirFlags );
    inputPathSubdirs( ismember( {inputPathSubdirs.name}, {'.', '..'} ) ) = [];

    % Check number of workers
    if options.workers <= 0
        error( 'aqua:incorrectWorkers', ...
               'Must have at least 1 worker for export (not %s)', options.workers );
    end

    % Close current parallel pool if it exists (for cleanup)
    poolobj = gcp( 'nocreate' );
    if numel( poolobj ) > 0
        delete( poolobj );
    end

    % Create parallel pool
    parpool( options.workers );

    % Iterate over each subdir
    parfor k = 1 : length( inputPathSubdirs )

        curDataset = inputPathSubdirs(k).name;

        disp( strcat( "Starting dataset #", num2str( k ), ": ", curDataset ) );

        curFilename = dir( fullfile( options.inputPath, curDataset, sprintf( "%s*.mat", curDataset ) ) ).name;
        curPath = fullfile( options.inputPath, curDataset, curFilename );

        try

            exportData = aqua_extract( curPath, 'metadata', metadata );

            saveFilename = strcat( curDataset, "_", options.saveSuffix, ".mat" );
            savePath = fullfile( options.outputPath, saveFilename );
            parsave( savePath, exportData );

            disp( strcat( "Completed dataset #", num2str( k ) ) );

        catch e

            disp( strcat( "Could not extract/save dataset #", num2str( k ) ) );
            disp( e );

        end

    end

    disp( 'Done!' );

end

%