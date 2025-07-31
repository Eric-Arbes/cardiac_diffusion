%This code provides a rough example on how to reconstruct obtained cardiac
%diffusion data. Please note that certain parameters like matrix size and
%partial Fourier factor are hardcoded in this example.
%Recon code is based on examples provided within Pulseq by Prof. Dr. Maxim
%Zaitsev.





%https://matlab.fandom.com/wiki/FAQ#How_can_I_process_a_sequence_of_files.3F
myFolder = 'C:\dummydir';  %location of rawdata files

seq_file_path = 'C:\dummydir\dummyseq.seq'; %Use ONE of the sequence files used to generated the data (i.e. just part 1 out of 3)

numberofaverages = 5;  %number of times any given sequence files was applied
numberofslices = 3;    %number of slices measured
numberofbreathholds = 3; %numbers of breathholds over which acquisition was split

slicestimesbreathholds = numberofslices*numberofbreathholds;
newmyFolderdel = append(myFolder,"\","temp","\*");   %assumes rawdata directory has nested subfolders "temp" and temp2" but this can be changed
newestmyFolderdel = append(myFolder,"\","temp\temp2","\*");
delete(newmyFolderdel)
delete(newestmyFolderdel)

filePattern = fullfile(myFolder, '*.dat');
theFiles = dir(filePattern);

%based on https://matlab.fandom.com/wiki/FAQ#How_can_I_process_a_sequence_of_files.3F
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fileindex = fullFileName(end-4);
    fileindex = str2num(fileindex);
    for o = 1:3

        data_file_path = fullFileName;
        
        slicenum = o;
        
        numim = 6; %number of separate kspaces contained in the data
        
        
        %run grappa sequence generation code to get linesuntilzerolines, NeoACS etc..
        Nyabs = 130;
        Nxref = 130;
        numim_orig = numim;
        numimabs = numim;
        twix_obj = mapVBVD(data_file_path);
        
        
        seq = mr.Sequence(); % Create a new sequence object
        seq.read(seq_file_path,'detectRFuse');
        [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
        
        if iscell(twix_obj)
            rawdata = double(twix_obj{end}.image.unsorted());
        else
            rawdata = double(twix_obj.image.unsorted());
        end
        rawdata_cal = rawdata(:,:,98*(slicenum-1)+1:98*(slicenum-1)+98);
        rawdata_ims = [];
        
        for i = slicenum:3:(numim*3)-3
            rawdata_ims_temp = rawdata(:,:,(98*3)+1+34*(i-1):(98*3)+34*(i));
            rawdata_ims = cat(3,rawdata_ims,rawdata_ims_temp);
        end
        
        rawdata = cat(3,rawdata_cal,rawdata_ims);
        
        numchan = size(rawdata,2);
        freqsize = size(rawdata,1);
        ACSbegin = 13;
        ACSend = 29;
        
        ACSbegin1 = 13;
        ACSend1 = 29;
        
        ACSbeginold = ACSbegin;
        ACSendold = ACSend;
        ACSdiff = ACSend-ACSbegin;
        ACSbegin = 25;
        ACSend = ACSbegin + ACSdiff;

        realindex = size(rawdata,3)- 98;
        realindex = realindex/(numim_orig-1);
        newrawdata = zeros(freqsize,numchan,numim*98); %hardcoded
        newrawdata(:,:,1:98) = rawdata(:,:,1:98);
        realis = 2:2:numim_orig;


        rawdatacopy = rawdata;
        

        %% automatic detection of the measurement parameters (FOV, matrix size, etc)
        nADC = size(rawdata, 1);
        k_last=ktraj_adc(:,end);
        k_2last=ktraj_adc(:,end-nADC);
        delta_ky=k_last(2)-k_2last(2);
        fov=1/abs(delta_ky);
        Ny_post=round(abs(k_last(2)/delta_ky));
        if k_last(2)>0
            Ny_pre=round(abs(min(ktraj_adc(2,:))/delta_ky));
        else
            Ny_pre=round(abs(max(ktraj_adc(2,:))/delta_ky));
        end
        Nx=2*max([Ny_post,Ny_pre]);
        Nx = Nxref;
        Ny_sampled=Ny_pre+Ny_post+1;
        
        
        rawdata = rawdata(:,:,1:98);
        rawcal = rawdata;
        
        
        %% classical phase correction / trajectory delay calculation 
        %  here we assume we are dealing with the calibration data
        data_odd=fftshift(ifft(fftshift(rawdata(:,:,1:2:end),1)),1);
        data_even=fftshift(ifft(fftshift(rawdata(end:-1:1,:,2:2:end),1)),1);
        
        if size(data_odd,3) > size(data_even,3)
            cmplx_diff=data_even.*conj(data_odd(:,:,1:size(data_even,3)));
        elseif size(data_odd,3) < size(data_even,3)
            cmplx_diff=data_even(:,:,1:size(data_odd,3)).*conj(data_odd);
        else
            cmplx_diff=data_even.*conj(data_odd);
        end
            
        cmplx_slope=cmplx_diff(2:end,:,:).*conj(cmplx_diff(1:end-1,:,:));
        mslope_phs=angle(sum(cmplx_slope(:)));
        dwell_time=(t_adc(nADC)-t_adc(1))/(nADC-1);
        measured_traj_delay=mslope_phs/2/2/pi*nADC*dwell_time;
        fprintf('measured trajectory delay (assuming it is a calibration data set) is %g s\n', measured_traj_delay);
        fprintf('type this value in the section above and re-run the script\n');

        traj_recon_delay = measured_traj_delay;
        rawdata = rawdatacopy;
        [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
        
        %% automatic detection of the measurement parameters (FOV, matrix size, etc)
        nADC = size(rawdata, 1);
        k_last=ktraj_adc(:,end);
        k_2last=ktraj_adc(:,end-nADC);
        delta_ky=k_last(2)-k_2last(2);
        fov=1/abs(delta_ky);
        Ny_post=round(abs(k_last(2)/delta_ky));
        if k_last(2)>0
            Ny_pre=round(abs(min(ktraj_adc(2,:))/delta_ky));
        else
            Ny_pre=round(abs(max(ktraj_adc(2,:))/delta_ky));
        end
        Nx=2*max([Ny_post,Ny_pre]);
        Nx = Nxref;
        Ny_sampled=Ny_pre+Ny_post+1;
        
        
        % %% classical phase correction / trajectory delay calculation 
        % %  here we assume we are dealing with the calibration data
        % % we do not calculate the constant phase term here because it depends on
        % % the definitions of the center of k-space and image-space 
        
        
        %% analyze the trajecotory, resample the data
        % here we expect rawdata ktraj_adc loaded (and having the same dimensions)
        nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
        nAcq = size(rawdata,3);
        nD=size(ktraj_adc, 1);
        
        kxmin=min(ktraj_adc(1,:));
        kxmax=max(ktraj_adc(1,:));
        kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
        kmaxabs=max(kxmax1, -kxmin);
        
        kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions
        ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
        t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);
        
        data_resampled=zeros(length(kxx), nCoils, nAcq);
        ktraj_resampled=zeros(nD, length(kxx), nAcq);
        t_adc_resampled=zeros(length(kxx), nAcq);
        
        for a=1:nAcq
            for c=1:nCoils
                data_resampled(:,c,a)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a),kxx,'spline',0);
            end
            ktraj_resampled(1,:,a)=kxx;
            for d=2:nD
                ktraj_resampled(d,:,a)=interp1(ktraj_adc2(1,:,a),ktraj_adc2(d,:,a),kxx,'linear',NaN);%used to be "linear"
            end
            t_adc_resampled(:,a)=interp1(ktraj_adc2(1,:,a),t_adc2(:,a),kxx,'linear',NaN);%used to be "linear"
        end
  
        
        %% in some cases because of the incorrectly calculated trajectory phase correction may be needed
        %  one such case is the use of the frequency shift proportional to gradient
        %  in combination with the gradient delay and FOV offset in the RO direction
        %  this calculation is best done with the calibration data, but also seems
        %  to work with the actual image data
        
        % here we assume we are dealing with the calibration data
        
        data_resampledcopy = data_resampled;
        data_resampled = data_resampled(:,:,1:98);
        
        data_odd=fftshift(ifft(fftshift(data_resampled(:,:,1:2:end),1)),1);
        data_even=fftshift(ifft(fftshift(data_resampled(:,:,2:2:end),1)),1);
        
        
        if size(data_odd,3) > size(data_even,3)
            cmplx_diff1=data_even.*conj(data_odd(:,:,1:size(data_even,3)));
        elseif size(data_odd,3) < size(data_even,3)
            cmplx_diff1=data_even(:,:,1:size(data_odd,3)).*conj(data_odd);
        else
            cmplx_diff1=data_even.*conj(data_odd);
        end
        if size(data_odd,3) > size(data_even,3)
            cmplx_diff2=data_even.*conj(data_odd(:,:,1:size(data_even,3)));
        elseif size(data_odd,3) < size(data_even,3)
            cmplx_diff2=data_even(:,:,1:size(data_odd,3)).*conj(data_odd);
        else
            cmplx_diff2=data_even.*conj(data_odd);
        end
        mphase1=angle(sum(cmplx_diff1(:)));
        mphase2=angle(sum(cmplx_diff2(:)));
        mphase=angle(sum([cmplx_diff1(:); cmplx_diff2(:)]));
        data_resampled = data_resampledcopy;
        %%
        pc_coef = mphase2/(2*pi);
        
        data_pc=data_resampled;
        for c=1:nCoils
            for i=1:size(data_resampled,1)
                data_pc(i,c,:)=squeeze(data_resampled(i,c,:)).*exp(1i*2*pi*pc_coef*((mod((1:size(data_pc,3))',2))));       
            end
        end
        
        %% reshape for multiple slices or repetitions
        
        %% display results
        
        
        % strictly speaking we have to do an intensity compensation here due to the
        % convolution at the interp1() step, but for now we ignore it...
        
        
        %%
        data_resampled = data_pc;
        
        raw = squeeze(data_resampled);
        
        realindex = size(raw,3)- 98;
        realindex = realindex/(numim_orig-1);
        newrawdata = zeros(130,numchan,numim*98); %hardcoded
        newrawdata(:,:,1:98) = raw(:,:,1:98);
        realis = 2:2:numim_orig;
        for i = 1:(numim-1)
            tempraw = raw(:,:,(99+(i-1)*realindex):98+i*realindex);
            tempnewraw = newrawdata(:,:,(99+(i-1)*98):98+i*98);
            if mod(AccelFac,2) == 0 && mod(i,2) == 0
            end
            counter = 1;
            if i == 1
                gigaindex = ceil((Ny/2)*partFourierFactor)-linesuntilzeroline;
                for j = gigaindex:gigaindex+NeoACS-1
                    tempnewraw(:,:,j) = tempraw(:,:,counter);
                    counter = counter + 1;
                end
            counter = 1;
            else    
                for j = 1:AccelFac:98
                    tempnewraw(:,:,j) = tempraw(:,:,counter);
                    counter = counter + 1;
                end
            end
            counter = 1;
            counter = 1;
            newrawdata(:,:,(99+(i-1)*98):98+i*98) = tempnewraw;
        end
        
        raw = newrawdata;
        
        raw = permute(raw,[1,3,2]);
        % only takes first reference scan at the moment
        % zero pad attempt
        padval = size(raw,1) - (size(raw,2)/numim);
        padmat = zeros(size(raw,1),1,padval);     
        rawdatatemp = [];
        rawdatabigtemp = [];
        for i = 1:size(raw,3)
            rawdatatemp = [];
            for j = 1:numim
                rawdatatemp_part = raw(:,1+(size(raw,2)/numim)*(j-1):(size(raw,2)/numim)*j,i);
                rawdatatemp_part = [squeeze(padmat) squeeze(rawdatatemp_part)];
                rawdatatemp = [rawdatatemp rawdatatemp_part];
            end
            rawdatabigtemp(:,:,i) = rawdatatemp;
        end
        
        raw = rawdatabigtemp; %zero filled!
        
        image = raw.*0;
        kspace = [];
        
        slicecolsize = size(raw,2)/numim;
        counter = 0;
        for i = 1:numim
            kspace_temp = zeros(size(raw,1),slicecolsize,size(raw,3));
            kspace_temp(:,1:slicecolsize,:) = raw(:,(counter*slicecolsize)+1:(counter*slicecolsize)+slicecolsize,:);
            kspace = cat(4,kspace,kspace_temp);
            counter = counter + 1;
        
        end
        kspace = permute(kspace,[1,2,4,3]);
        testspace = squeeze(kspace(:,:,2,:));
        testspace = permute(testspace,[3,2,1]);
        testspace = permute(testspace,[1,3,2]);
        
        testspaceunder = squeeze(kspace(:,:,3,:));
        testspaceunder = permute(testspaceunder,[3,2,1]);
        testspaceunder = permute(testspaceunder,[1,3,2]);
        
        Ny = Nyabs;
        testspaceref = testspace(:,:,gigaindex+(Ny-98):gigaindex+NeoACS-1+(Ny-98));
        R       =   [1,AccelFac];
        kernel  =   [3,4];
        
        recon   =   grappa(testspaceunder, testspaceref, R, kernel);
        show_quad(testspaceunder, recon,'grappa');
        
        newkspace = [];
        weights = [];
        counter = 1;
        
        %   grappa.m by Mark Chiew https://github.com/mchiew/grappa-tools
        %   mchiew@fmrib.ox.ac.uk
        
        
        for i = 1:numim-2
            testspace = squeeze(kspace(:,:,2,:));
            testspace = permute(testspace,[3,2,1]);
            testspace = permute(testspace,[1,3,2]);
            
            testspaceunder = squeeze(kspace(:,:,i+2,:));
            testspaceunder = permute(testspaceunder,[3,2,1]);
            testspaceunder = permute(testspaceunder,[1,3,2]);
            Ny = Nyabs;
            testspaceref = testspace(:,:,gigaindex+(Ny-98):gigaindex+NeoACS-1+(Ny-98));
            R       =   [1,AccelFac];
            kernel  =   [3,4];
            
            recon   =   grappa(testspaceunder, testspaceref, R, kernel);
            weight = grappaweights(testspaceunder, testspaceref, R, kernel);
            newkspace(counter,:,:,:) = recon;
            weights(counter,:,:) = weight;
            counter = counter+1;
        end
        
        %from Dr. Johannes Fischer
        windowLength = 6;
        hannWindow = hann(windowLength);
        halfwindow = windowLength/2;
        
        newkspace = permute(newkspace,[3,4,1,2]);
        %from Dr. Johannes Fischer
        image = fft2stack(kspace);
        
        newimage = fft2stack(newkspace);
        image_sos = abs(sqrt(sum(abs(image).^2,4)));
        newimage_sos = abs(sqrt(sum(abs(newimage).^2,4)));
        meanim = image_sos(:,:,2:end);
        meanim = mean(meanim,3);



    str = [myFolder,"\",num2str(k),num2str(fileindex),num2str(o),".mat"];
    str = append(myFolder,"\","temp","\",num2str(k),num2str(fileindex),num2str(o),".mat");
    save(str,"newimage_sos");
    end
end

close all


newmyFolder = append(myFolder,"\","temp");
newfilePattern = fullfile(newmyFolder, '*.mat');
newtheFiles = dir(newfilePattern);

newFileNames = {newtheFiles.name};
newFileNumbers = [];
for i = 1:length(newFileNames)
    tempnewname = newFileNames{i};
    tempnewname = strsplit(tempnewname,'.');
    tempnewname = tempnewname{1};
    tempnewnumbers = str2num(tempnewname);
    newFileNumbers(end+1) = tempnewnumbers;
end

 

%https://de.mathworks.com/matlabcentral/answers/67235-sorting-the-name-field-in-dir-command
[~,ind] = sort(newFileNumbers);

sortednewtheFiles = newtheFiles(ind);
newtheFiles = sortednewtheFiles;

%https://matlab.fandom.com/wiki/FAQ#How_can_I_process_a_sequence_of_files.3F
for j = 1:slicestimesbreathholds
    fullfiles = [];
    for k = j:slicestimesbreathholds:length(newtheFiles)
        newbaseFileName = newtheFiles(k).name;
        newfullFileName = fullfile(newtheFiles(k).folder, newbaseFileName);
        newfullFile = load(newfullFileName);
        newfullFile = newfullFile.newimage_sos; 
        fullfiles = cat(4,fullfiles,newfullFile);
    end
    averaged = sum(fullfiles,4);
    averaged = averaged/numberofaverages;
    str = append(myFolder,"\","temp","\","temp2","\",num2str(j),".mat");
    save(str,"averaged");
end

newestmyFolder = append(newmyFolder,"\","temp2");
newestfilePattern = fullfile(newestmyFolder, '*.mat');
newesttheFiles = dir(newestfilePattern);

dellist = [];

%removes unwanted files, in this example case tSNR data that was saved to
%the same directory
for i = 1:length(newesttheFiles)
    tempsortname = newesttheFiles(i).name;
    if contains(tempsortname,'tSNR')
        dellist(end+1) = i;
    end
end
newesttheFiles(dellist) = [];





newestFileNames = {newesttheFiles.name};
newestFileNumbers = [];
for i = 1:length(newestFileNames)
    tempnewname = newestFileNames{i};
    tempnewname = strsplit(tempnewname,'.');
    tempnewname = tempnewname{1};
    tempnewnumbers = str2num(tempnewname);
    newestFileNumbers(end+1) = tempnewnumbers;
end

 

%https://de.mathworks.com/matlabcentral/answers/67235-sorting-the-name-field-in-dir-command
[~,ind] = sort(newestFileNumbers);

sortednewesttheFiles = newesttheFiles(ind);
newesttheFiles = sortednewesttheFiles;

finalslices = {};
for k = 1:numberofslices
    newestfullfiles = {};
    for i = k:numberofbreathholds:slicestimesbreathholds
            newestbaseFileName = newesttheFiles(i).name;
            newestfullFileName = fullfile(newesttheFiles(i).folder, newestbaseFileName);
            newestfullFile = load(newestfullFileName);
            newestfullFile = newestfullFile.averaged;
            newestfullfiles{end+1} = newestfullFile;
    end
    %hardcoded for now assuming 3 slices
    tempslice1 = newestfullfiles{1};
    tempslice2 = newestfullfiles{2};
    tempslice3 = newestfullfiles{3};
    slicex = cat(3,tempslice1,tempslice2(:,:,2:end));
    slicex = cat(3,slicex,tempslice3(:,:,2:end));
    finalslices{end+1} = slicex;
end

slice1 = finalslices{1};
slice2 = finalslices{2};
slice3 = finalslices{3};

%%
save("C:\dummydir\slice1_example.mat","slice1");
save("C:\dummydir\slice2_example.mat","slice2");
save("C:\dummydir\slice3_example.mat","slice3");


