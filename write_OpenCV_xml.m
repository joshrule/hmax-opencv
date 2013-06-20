function write_OpenCV_xml(filename,struct)

%History: 07/12/12
%Code to write xml code that imitates what OpenCV XML/YAML persistence code
%(writer only). My code is similar to existing write_UAV_par.m for
%another program.
%For now assume this "double" notation: 1.0034550000000000e+02
%Updates:
%07/13/2012: It can handle scalars, vectors and multi-dimensional matrix
%07/15/2012: It recognizes if an element in the structure is a string and
%it will written as such.

%Typical Inputs:
%1. filaname = should end with ".xml"
%2. structure = can only contain matrices (of any dimension)
%                and be one layer deep down.

%Typical Input:
% struct.a=magic(5);
% struct.b=magic(6);
% struct.c=[1,2,3,4;5,6,7,8];
% struct.d=rand(2,3,4,5);
%write_OpenCV_xml('test_matlab_3.xml',struct)


%Initialization
format='%1.16e\t';%'%1.16e\t';  %May make this an input
line_ele=3; %number of elements per line


%identify data in structure:
names = fieldnames(struct);
num_files=length(names); %really number of variables

%Open file
hpar=fopen(filename,'w');

%Initial header:
fprintf(hpar, '%s\n', '<?xml version="1.0"?>');
fprintf(hpar, '%s\n', '<opencv_storage>');


for ij=1:num_files
    
    %Read data
    value = getfield(struct, names{ij});
    
    if ischar(value) %07/15/2012
      fprintf(hpar, '%s\n', strcat('<',names{ij},'>',value,'</',names{ij},'>'));

        
    else
        
        sz=size(value);
        
        
        %Rearrange data from column major to row major
        if length(sz) > 2
            %         new_value=value.*0;
            %         num_2d_matrix=prod(sz(3:end));
            %
            %         for jh=1:num_2d_matrix
            %             new_value(:,:,jh)=value(:,:,jh).';
            %         end
            
            %OR using Matlab function:
            new_value=permute(value,[2 1 3:length(sz)]);
            
            value_t=new_value;
            
        else
            value_t=value.';
        end
        
        
        %Choose field information based on original dimensionality of data:
        if length(sz) == 2
            fprintf(hpar, '%s\n', strcat('<',names{ij},' type_id="opencv-matrix">'));
            fprintf(hpar, '%s\n', strcat('  <rows>',num2str(sz(1)),'</rows>'));
            fprintf(hpar, '%s\n', strcat('  <cols>',num2str(sz(2)),'</cols>'));
        elseif length(sz) == 3
            fprintf(hpar, '%s\n', strcat('<',names{ij},' type_id="opencv-nd-matrix">'));
            fprintf(hpar, '%s\n', ['  <sizes>' num2str(sz(1)) ' ' num2str(sz(2)) ' ' num2str(sz(3)) '</sizes>']);
        else %dimension higher than 3
            fprintf(hpar, '%s\n', strcat('<',names{ij},' type_id="opencv-nd-matrix">'));
            fprintf(hpar, '%s','  <sizes>');
            for jh=1:length(sz)
                fprintf(hpar, '%s',[num2str(sz(jh)) ' ']);
            end
            fprintf(hpar, '%s\n','</sizes>');
            
        end
        fprintf(hpar, '%s\n', strcat('  <dt>','d','</dt>'));
        fprintf(hpar, '%s\n', '  <data>');
        
        
        
        
        %1. One way:
        %Put for loop here and write 3 numbers at the time, followed by \r or
        %\r\n
        ser_data=value_t(:);
        num_ele=length(ser_data);
        num_loop=ceil(num_ele/line_ele);
        
        ct1=1;
        for kl=1:num_loop
            if kl == num_loop
                a_str = sprintf(format,ser_data(ct1:end));%or value(:)
            else
                a_str = sprintf(format,ser_data(ct1:ct1+line_ele-1));%or value(:)
            end
            if ispc
                a_str = strrep(a_str, 'e+0', 'e+');
            end
            ct1=ct1+line_ele;
            fprintf(hpar, '%s\n', a_str);
        end
        
        %     %2. Second way
        %     %Change data to strings
        %     a_str = sprintf(format,value_t(:));%or value(:)
        %     %Write serialized data
        %     fprintf(hpar, '%s\n', a_str);
        
        %finalize matrix:
        fprintf(hpar, '%s\n', '  </data>');
        fprintf(hpar, '%s\n', strcat('</',names{ij},'>'));
        
    end
    
end

%Final header:
fprintf(hpar, '%s\n', '</opencv_storage>');

%Close file
fclose(hpar);



end
