function parsave( fn, saveStruct )
%PARSAVE Save `saveStruct` to path `fn` as an HDF5-compatible file
%   (Implemented so as to be parallelizeable.)

save( fn, '-struct', 'saveStruct', '-v7.3' );

end

%