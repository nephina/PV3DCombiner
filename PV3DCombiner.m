clear
clc


FileSearch = dir('PV3D Files');
FileTable = struct2cell(FileSearch);

FileNameList = FileTable(1,3:end);
NumberofFiles = length(FileNameList);


for i = 1:NumberofFiles
    file = FileNameList{i};
    FileName = sprintf('%c',file(1,:));
    FileNames(i,:) = sprintf(FileName);
end

File = zeros(1,10);
%% Read .csv data from each file
for i = 1:NumberofFiles
        loadFile =  sprintf('PV3D Files/%s', FileNames(i,:));
        tempfile = dlmread(loadFile,",", 2, 0);
        File = cat(1,File,tempfile);
end
File = File(2:end,:);

X = File(:,1); Y = File(:,2); Z = File(:,3); U = File(:,4); V = File(:,5); W = File(:,6);
Ustd = std(U); Vstd = std(V); Wstd = std(W);
Umean = mean(U); Vmean = mean(V); Wmean = mean(W);

VelocityAbs = sqrt(((File(:,4)).^2)+((File(:,5)).^2)+((File(:,6)).^2));

figure(1,'position',get(0,'screensize'));
hold on
subplot(2,3,1);
hist(File(:,4),1000);
%set(gca,'xscale','log');
axis('auto');
subplot(2,3,2);
hist(File(:,5),1000);
%set(gca,'xscale','log');
axis('auto');
subplot(2,3,3);
hist(File(:,6),1000);
%set(gca,'xscale','log');
axis('auto');
subplot(2,3,4:6);
hist(VelocityAbs,1000);
%set(gca,'xscale','log');
axis('auto');

FilterParameter = input('What is the standard deviation filter level?: '); % 0.01*input('What is the filter level? (% of velocity range below which particles are kept)');



VelocityMask = ones(size(File,1),1);
VelocityMask(U > ((FilterParameter*Ustd)) | U < (-(FilterParameter*Ustd)) | V > ((FilterParameter*Vstd)) | V < (-(FilterParameter*Vstd)) | W > ((FilterParameter*Wstd)) | W < (-(FilterParameter*Wstd))) = 0;
NumofFilteredPoints = sum(VelocityMask);


%% Kernel Filtering
KernelSize = 0.01*input('What is the kernel size as percentage of the whole?: ');
KernelBounds = [min(X):(0.25*range(X)*KernelSize):max(X);min(Y):(0.25*range(Y)*KernelSize):max(Y);min(Z):(0.25*range(Z)*KernelSize):max(Z)];
Points = [X Y Z];
Velocities = [U V W];
OtherData = [File(:,7) File(:,8) File(:,9)];

FilteredFile = [];
for i = 1:(4/KernelSize)
    for j = 1:(4/KernelSize)
        for k = 1:(4/KernelSize)
            Kernelindices = Points(:,1) > KernelBounds(1,i) & Points(:,1) < KernelBounds(1,i+1) & Points(:,2) > KernelBounds(2,j) & Points(:,2) < KernelBounds(2,j+1) & Points(:,3) > KernelBounds(3,k) & Points(:,3) < KernelBounds(3,k+1);
            KernelPoints = Points(Kernelindices,:);
            if length(KernelPoints) > 2
                KernelVelocities = Velocities(Kernelindices,:);
                KernelOtherData = OtherData(Kernelindices,:);
                KernelAll = cat(2,KernelPoints,KernelVelocities,KernelOtherData);
                KernelStd = std(KernelVelocities,0,1);
                KernelMean = mean(KernelVelocities,1);
                KernelMedian = median(KernelVelocities,1);
                KernelMAD = 2*median(abs(KernelVelocities-KernelMedian));
                if KernelMAD == 0
                else
                    KernelFilter = KernelVelocities(:,1) < (KernelMedian(1)+(KernelMAD(1)*FilterParameter)) & KernelVelocities(:,1) > (KernelMedian(1)-(KernelMAD(1)*FilterParameter)) & KernelVelocities(:,2) < (KernelMedian(2)+(KernelMAD(2)*FilterParameter)) & KernelVelocities(:,2) > (KernelMedian(2)-(KernelMAD(2)*FilterParameter)) & KernelVelocities(:,3) < (KernelMedian(3)+(KernelMAD(3)*FilterParameter)) & KernelVelocities(:,3) > (KernelMedian(3)-(KernelMAD(3)*FilterParameter));
                    FilteredFile = cat(1,FilteredFile,KernelAll(KernelFilter,:));
                end
            else
            end
        end
    end
end

NumofFilteredPoints = length(FilteredFile);

filetype = menu('Which file type do you want?','Combined PV3D','.xyz pointcloud','.ply vertex and normals','exit');
switch filetype
    case 1 % Write .pv3d file for V3V processing
        filename = fopen('Combined.pv3d','w');
        header = sprintf('Title="D:\ExperimentsV3V\MCW09Alexa\May7th\Analysis\bgrem_May7th000001.T000.D000.P000.H000.pv3d" VARIABLES="X","Y","Z","U","V","W","CHC","idParticleMatchA","idParticleMatchB",DATASETAUXDATA DataType="P",DATASETAUXDATA Dimension="3",DATASETAUXDATA HasVelocity="Y",DATASETAUXDATA ExtraDataNumber="2",ZONE T="T1",I=%d,F=POINT,\n',length(FilteredFile));
        fprintf(filename,header);
        for i = 1:length(FilteredFile)
            fprintf(filename,'\n%d,%d,%d,%d,%d,%d,%d,%d,%d,',FilteredFile(i,1),FilteredFile(i,2),FilteredFile(i,3),FilteredFile(i,4),FilteredFile(i,5),FilteredFile(i,6),FilteredFile(i,7),FilteredFile(i,8),FilteredFile(i,9))
        end
        fclose(filename);
    case 2 % Write .xyz file
        filename = fopen('UnstructPointCloud.xyz','w');
        for i = 1:length(FilteredFile)
            fprintf(filename,'%d %d %d\n',FilteredFile(i,1),FilteredFile(i,2),FilteredFile(i,3));
        end
        fclose(filename);
    case 3 % Write .ply file
        filename = fopen('UnstructPointNormals.ply','w');
        fprintf(filename,'ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\nend_header\n',NumofFilteredPoints);

        for i = 1:length(FilteredFile)
            fprintf(filename, '%d %d %d %d %d %d\n',FilteredFile(i,1),FilteredFile(i,2),FilteredFile(i,3),FilteredFile(i,4),FilteredFile(i,5),FilteredFile(i,6));
            %if VelocityMask(i) ~= 0       % VelocityAbs(i) < FilterParameter*MaxVelocity
            %    fprintf(filename, '%d %d %d %d %d %d\n',File(i,1),File(i,2),File(i,3),File(i,4),File(i,5),File(i,6));
            %end
        end
        fclose(filename);
    case 4
end