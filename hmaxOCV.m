function hmaxOCV(imgList,patchFile,maxSize,nProcs)
% hmaxOCV(imgList,patchFile,maxSize,nProcs)
%
% generate a set of HMAX-OCV activations
%
% imgList: a cell array, the names of images to process
% patchFile: a char array (string), the name of the XML file holding patch info
% nProcs: the number of processes to use (i.e. number of processors)
%
% note: requires 'cp','find','xargs','mkdir', and 'echo'

    if (nargin < 4) [~,sysResult] = system('nproc'); nProcs = str2num(sysResult); end;

    % make nProcs directories
    system('rm -rf /home/joshrule/maxlab/opencv/images/*');
    dirNames = arrayfun(@(x) ['/home/joshrule/maxlab/opencv/images/' num2str(x) '/'],...
                        0:nProcs-1,'UniformOutput',0);
    cellfun(@(x) system(['mkdir ' x]),dirNames,'UniformOutput',0);
    dirString = [sprintf('%s\\n',dirNames{1:end-1}), dirNames{end}];

    % split all images across nProcs directories
    for iImg = 1:length(imgList)
        imgDirs(iImg) = mod(iImg,nProcs);
        [fa,fb,fc] = fileparts(imgList{iImg});
        newImgFile{iImg} = ['/home/joshrule/maxlab/opencv/images/' ...
                      num2str(imgDirs(iImg)) '/' fb lower(fc)];
        while ~exist(newImgFile{iImg},'file')
            imwrite(uint8(resizeImage(double(imread(imgList{iImg})), ...
                                      maxSize)), newImgFile{iImg}, 'Quality',100);
        end
    end

    % run the simulation
    [pa,pb,pc] = fileparts(patchFile);
    system(['echo -e "' dirString '" | xargs -I % -n 1 -P 0'...
            ' ~/maxlab/opencv/bin/HMAX_FE -id % -pd ' pa '/ -pf ' pb pc...
            ' > ~/maxlab/opencv/hmaxOCV.log']);

    % copy the results back!
    for iImg = 1:length(imgList)
        [fa,fb,fc] = fileparts(imgList{iImg});
        system(['cp /home/joshrule/maxlab/opencv/images/'...
                num2str(imgDirs(iImg)) '/' fb '.xml ' fa]);
    end
end
