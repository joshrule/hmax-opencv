function xmlFile = matlabPatches2OCVPatches(filters,filterSizes,c1Scale, ...
  c1Space,c1OL,patches,patchSizes,patchName,patchDir)
% xmlFile = matlabPatches2OCVPatches(filters,filterSizes,c1Scale, ...
%  c1Space,c1OL,patches,patchSizes,patchName,patchDir)
    matFile = [patchDir patchName '.mat'];
    xmlFile = [patchDir patchName '.xml'];
    s = matlabPatches2aplPatches(filters,filterSizes,c1Scale,c1Space,c1OL, ...
      patches,patchSizes);
    truePatchName = genvarname(patchName);
    dummy.(truePatchName) = s;
    save(matFile,'-struct','dummy');
    aplPatches2xmlPatches(xmlFile,matFile,truePatchName);
end 
