function aplPatches2xmlPatches(xmlFile,matFile,patchSet)
% aplPatches2xmlPatches(xmlFile,matFile,patchSet)
%
% convert APL patch struct to XML patch file
%
% xmlFile: the resultant XML file
% matFile: the file holding the APL patch struct
% patchSet: the name of the APL patch struct itself.
    matVars = load(matFile); % load Mat file

    patches = getfield(matVars, patchSet); %verify object name
    patches.suppressionparams = patches.s1c1suppress;
    patches.pruningparams = 0;
    patches.kernel = 'linear';

    write_OpenCV_xml(xmlFile,patches);
end
