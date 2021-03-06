function hmaxOCV(imgList,patchFile,hmaxHome,maxSize,nProcs)
% hmaxOCV(imgList,patchFile,hmaxHome,maxSize,nProcs)
%
% generate a set of HMAX-OCV activations
%
% imgList: a cell array, the names of images to process
% patchFile: a char array (string), the name of the XML file holding patch info
% hmaxHome: string, top-level directory of HMAX-OCV repository
% nProcs: the number of processes to use (i.e. number of processors)
%
% note: requires 'cp','find','xargs','mkdir', and 'echo'

    if (nargin < 5) [~,result] = system('nproc'); nProcs = str2num(result); end;

    % make nProcs directories
    system(['rm -rf ' hmaxHome 'images/*']);
    dirNames = arrayfun(@(x) [hmaxHome 'images/' num2str(x) '/'],...
                        0:nProcs-1,'UniformOutput',0);
    cellfun(@(x) system(['mkdir ' x]),dirNames,'UniformOutput',0);
    dirString = [sprintf('%s\\n',dirNames{1:end-1}), dirNames{end}];

    % split all images across nProcs directories
    for iImg = 1:length(imgList)
        imgDirs(iImg) = mod(iImg,nProcs);
        [fa,fb,fc] = fileparts(imgList{iImg});
        newImgFile{iImg} = [hmaxHome 'images/' num2str(imgDirs(iImg)) '/' ...
          fb '.bmp'];
        while ~exist(newImgFile{iImg},'file')
            system(['cp ' imgList{iImg} ' ' newImgFile{iImg}]);
            try 
                imwrite(uint8(resizeImage(double(imread(imgList{iImg})),maxSize)), ...
                    newImgFile{iImg},'bmp');
            catch
                system(['convert ' imgList{iImg} ' -colorspace rgb ' imgList{iImg}]);
                imwrite(uint8(resizeImage(double(imread(imgList{iImg})),maxSize)), ...
                    newImgFile{iImg},'bmp');
            end
        end
    end

    % run the simulation
    [pa,pb,pc] = fileparts(patchFile);
    system(['echo -e "' dirString '" | xargs -I % -n 1 -P 0 ' hmaxHome ...
      'bin/HMAX_FE -id % -pd ' pa '/ -pf ' pb pc ' > ' hmaxHome 'hmaxOCV.log']);

    % copy the results back!
    for iImg = 1:length(imgList)
        [fa,fb,fc] = fileparts(imgList{iImg});
        system(['cp ' hmaxHome 'images/' num2str(imgDirs(iImg)) ...
          '/' fb '.xml ' fa]);
    end
end
