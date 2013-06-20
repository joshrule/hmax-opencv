function [out_struct]=read_OpenCV_xml(filename,varargin)
%Reading the xml files, I just created with write_OpenCV_xml.m
%Created by Pedro A. Rodriguez

if ~isempty(varargin)
   ignoreText = varargin{1};
else
   ignoreText = 1;
end

%Typical input:
% struct=read_OpenCV_xml('classifierstatefile_xml_opencv_svmStruct.xml')
out_struct=[];
% xDoc = xmlread(filename);
struct = parseXML(filename);

num_fields=length(struct.Children);
%idx=2:2:num_fields;

cnt = 0;
for ij=1:length(struct.Children)
    if ignoreText && (isequal(struct.Children(ij).Name,'#comment') || ...
        isequal(struct.Children(ij).Name,'#text'))
        continue;
    else
        cnt = cnt + 1;       
    end

    out_struct(cnt).name=struct.Children(ij).Name;
    if ~isempty(struct.Children(ij).Attributes)
        out_struct(cnt).type=struct.Children(ij).Attributes.Value;
    else
        %assume string
        out_struct(cnt).type='string';
    end
    
    child=struct.Children(ij).Children;
    
    switch out_struct(cnt).type
        
        case 'opencv-matrix'
            n_rows=str2num(child(2).Children.Data);
            n_cols=str2num(child(4).Children.Data);
            data=child(8).Children.Data; %data in string
            if ispc
                data = strrep(data, 'e+','e+0');
            end
            dt=child(6).Children.Data;
            
            if length(dt) == 1
                %change data to matrix
                %             A=sscanf(data,['%' dt],[n_rows,n_cols]);
                %             A=sscanf(data,'%f',[n_rows,n_cols]);
                %             out_struct(ij).data=A;
                
                %Maybe I need to transpose it:
                A=sscanf(data,'%f',[n_cols,n_rows]);
                out_struct(cnt).data=A.';
                
            else
                %Create special case if multiple channels:
                num_channels=str2num(dt(2:end-2));
                A=sscanf(data,'%f',n_cols*n_rows*num_channels);
                A_reshape=reshape(A,[num_channels, n_cols n_rows]);
                out_struct(cnt).data=permute(A_reshape,[3 2 1]);
                
                
            end
        case 'string'
            if isequal(struct.Children(ij).Name,'#comment') || ...
                    isequal(struct.Children(ij).Name,'#text')
                out_struct(cnt).data=struct.Children(ij).Data;
            else
                out_struct(cnt).data=child.Data;
            end
            
        case 'opencv-nd-matrix'
            sizes=str2num(child(2).Children.Data);
            dt=child(4).Children.Data;
            data=child(6).Children.Data; %data in string
            if ispc
                data = strrep(data, 'e+','e+0');
            end
            
            %Maybe I need to transpose it:
            A=sscanf(data,'%f',prod(sizes));
            A_reshape=reshape(A,[sizes(2) sizes(1) sizes(3:end)]);
            out_struct(cnt).data=permute(A_reshape,[2 1 3:length(sizes)]);
            
        otherwise
            error('Failed to recognize XML type %s.',filename);
            
            
    end
    
end

end


function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems
% with very deeply nested trees.
try
    theStruct = parseChildNodes(tree);
catch
    error('Unable to parse XML file %s.',filename);
end
end

% ----- Subfunction PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1, numChildNodes);
    
    children = struct(             ...
        'Name', allocCell, 'Attributes', allocCell,    ...
        'Data', allocCell, 'Children', allocCell);
    
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end
end

% ----- Subfunction MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
    'Name', char(theNode.getNodeName),       ...
    'Attributes', parseAttributes(theNode),  ...
    'Data', '',                              ...
    'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end
end

% ----- Subfunction PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end
end