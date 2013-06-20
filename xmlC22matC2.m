function c2 = xmlC22matC2(imgNames,patchSet)
% c2 = xmlC22matC2(imgNames,patchSet)
%
% read C2 responses for image files
%
% imgNames: cell array, names of images
% patchSet: string, the patchSet for which to read activations
%
% c2: array, a standard C2 matrix
	activations = cell(length(imgNames), 1);
	parfor iImg = 1:length(imgNames)
        [pathstr, name, ~] = fileparts(imgNames{iImg});
        xmlName = fullfile(pathstr,[name '.xml']);
        tagMatch = findTagsInXMLFile(xmlName,{patchSet});
        activations{iImg} = tagMatch.(patchSet);
	end
    c2 = cell2mat(activations)';
end

function tagMatches = findTagsInXMLFile(xmlName,tags)
% classifierstate = findTags(xmlName,tags)
%
% Reads specified tags in a single xml file.
%
% xmlName: string, the name of the XML file
% tags: cell array of strings, the tags to read from the file
%
% tagMatches: the read tags
    fileStruct=read_OpenCV_xml(xmlName);

    nFields=length(fileStruct);
    nTags=length(tags);

    for iTag=1:nTags
        for iField=1:nFields
            nami=fileStruct(iField).name;
            out_data=fileStruct(iField).data;
            
            if strcmp(nami,tags{iTag})
                tagMatches.(tags{iTag}) = out_data;
            end
        end
    end
end

