function s = matlabPatches2aplPatches(filters,filterSizes,c1Scale,c1Space,c1OL,patches,patchSizes)
% s = matlabPatches2aplPatches(filters,filterSizes,c1Scale,c1Space,c1OL,patches,patchSizes)
%
% convert HMAX-MATLAB patch information to an APL patch struct
%
% filters, filterSizes, c1Space, c1OL: see HMAX-MATLAB's "C1.m"
% patches: cell array, each patches{i} is an array of patches as used by 
%     HMAX-MATLAB's "C2.m"
% patchSizes: array, patchSizes(1,i) is the number of rows in each orientation
%     of a patch in patches{i}. patchSizes(2,i) is the number of columns.
%
% s: struct, an APL patch struct
    s.filters = filters;
    s.fSiz = filterSizes;
    s.c1ScaleSS = c1Scale;
    s.c1SpaceSS = c1Space;
    s.c1OL = c1OL;
    s.s1c1suppress = 0;
    s.num_patches = length(patches);
    for iSize = 1:s.num_patches
        s = setfield(s,['cPatches_' num2str(iSize)],patches{iSize});
    end
    s.patchSizes = patchSizes(1:2,:);
end
