clear all;
close all;
tic

%SCRIPT: dataset: T2_star (Template customized for 10 TEs; 128x128 inplane X&Y FOV)
%LAST MODIFIED: 12/15/17 Adisetiyo
%        convert dicoms to nifti format
%        coregister individual TEs to first TE
%        create noise roi nifti file
%        calculate R2* fit to create map (1/sec)

list_subjects = {'S5a_post'}; % list of subjects IDs, ex: {'SubjectID' 'SubjectID'};
folder_base = 'P:\Play_Area\vadisetiyo\Projects\KTGF_Project\Cases\Data_MRI\Processed'; % folder where all case folders are located
spm8_path = 'H:\Programs\spm8_r4667'; 
addpath(spm8_path);


for a=1:length(list_subjects);

%--------------------------------------------------------------------------
% SPM8: convert dicom images to nifti format  (NOTE: from MFC_preprocess)
%--------------------------------------------------------------------------
  
 mkdir(fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\SPM_nifti'));        
 folder_dcm=fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\TE_all'); % work with either one input dcm folder or separated TE folders (just need make loop)
 folder_nii=fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\SPM_nifti');
   
    % get list of images from folder
    P = spm_select('fplist', folder_dcm, '.*');  % fplist: create list w/ full path to each dcm file within dcm_folder

    % read dicom headers
    hdr = spm_dicom_headers(P);

    % convert dicom to nii (SPM8 defualt naming)
    current_path = pwd;
    cd(folder_nii)
    fn_list_out = spm_dicom_convert(hdr, 'all', 'flat', 'nii'); 
    fn_list_out = fn_list_out.files; % path list of all output .nii files
    cd(current_path)
    
    % rename nii files accroding to TE as reflected in original .nii filename 
    for i = 1:length(fn_list_out) 
     if ~isempty(fn_list_out{i})
          idx1 = strfind(fn_list_out{i}, '-'); % index number of '-' in string of full file path to each .nii file 
          idx2 = strfind(fn_list_out{i}, '.');
          I(i) = str2num(fn_list_out{i}((idx1(end)+1):(idx2(end)-1)));  % take the number between the last '-' and the last '.'  
          if (I(i))<10
              fn_new = fullfile(folder_nii, ['TE_0' num2str(I(i)) '.nii']);
          else
              fn_new = fullfile(folder_nii, ['TE_' num2str(I(i)) '.nii']);
          end
          movefile(fn_list_out{i}, fn_new)
     else
            error('NIfTI files were not created! Perhaps path to DICOM images was incorrect.')
     end
    end
    
%--------------------------------------------------------------------------
% co-register all images to first TE magnitude images  (NOTE: from MFC_preprocess)
%--------------------------------------------------------------------------

  idx_TEs={'01' '02' '03' '04' '05' '06' '07' '08' '09' '10'};
  fn_target = fullfile(folder_nii,'TE_01.nii'); % file to register all images to
    
    for b=1:length(idx_TEs); % coregester each TE to the first TE one by one in loop
    fn_source=fullfile(folder_nii,['TE_' , idx_TEs{b}, '.nii']); % file to register to target in order to get registration matrix
    fn_other=fullfile(folder_nii,['TE_' , idx_TEs{b}, '.nii']);  % file to get registration matrix applied to
        
        % coregistration and reslicing parameters
        estflg.cost_fun = 'nmi';
        estflg.sep      = [4 2];
        estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        estflg.fwhm     = [7 7];
        wrtflg.interp   = 1;
        wrtflg.wrap     = [0 0 0];
        wrtflg.mask     = 0;

        hdr_trg = spm_vol(fn_target);
        hdr_src = spm_vol(fn_source);

        x  = spm_coreg(hdr_trg, hdr_src, estflg);
        M  = inv(spm_matrix(x));
        MM = zeros(4, 4, size(fn_other, 1));
        for j = 1:size(fn_other, 1)
            MM(:,:,j) = spm_get_space(fn_other(j,:));
        end
        for j = 1:size(fn_other, 1)
            spm_get_space(fn_other(j,:), M*MM(:,:,j));
        end

        fn_other = str2mat(fn_target, fn_other);
        spm_reslice(fn_other, wrtflg);

    end

%--------------------------------------------------------------------------
% automatically create noise ROI file 
%--------------------------------------------------------------------------
    
    file_TE1 = fullfile(folder_nii,'TE_01.nii');
    cd(folder_nii);
    
    HeaderInfo_TE1 = spm_vol(file_TE1);
    Data_TE = spm_read_vols(HeaderInfo_TE1);
    
    Data_TE1 = zeros(size(Data_TE));
    Data_TE1(:,:,:) = 0;   % avoid phase encoding direction, use TE_01 hdr info
    Data_TE1(4:11,4:124,:)= 1; % customized for 128x128 in plane FOV
    Data_TE1(118:125,4:124,:)= 1; % customized for 128x128 in plane FOV
   
    HeaderInfo_TE1.fname = [list_subjects{a},'_noise.nii']; 
    HeaderInfo_TE1.private.dat.fname = HeaderInfo_TE1.fname;
    spm_write_vol(HeaderInfo_TE1,Data_TE1);
     
end

for a=1:length(list_subjects);
       
%--------------------------------------------------------------------------
% calculate R2 fits 
%--------------------------------------------------------------------------

    clock
    tic
    
    TE=[4.92 9.84 14.76 19.68 24.6 29.52 34.44 39.36 44.28 49.20]';

    mkdir(fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\R2star_map_nlin'));
    folder_maps = fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\R2star_map_nlin');
    folder_noise = fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\SPM_nifti',[list_subjects{a}, '_noise.nii']);

    % make list of datasets
    cd(fullfile(folder_base,[list_subjects{a}],'dicom\T2_star\SPM_nifti'));
    list_files=dir('r*.nii'); % all input coregistered .nii files

    set_num=size(list_files,1); 
    F=strvcat(list_files.name);

    % initialize a header by loading an existing one and editing it
    temp=spm_vol(F(1,:)); %1st nifti
    hdr_temp=temp;

    % load images in as 4D matrix
    dat_echo=spm_vol(F);  % load header info from all TE nifti files 

    xdim=dat_echo(1).dim(1); % refer to first file in 5X1 structure array 'echo_dat' then refer to first component of dimenion field 'dim'
    ydim=dat_echo(1).dim(2);
    slnum=dat_echo(1).dim(3);
 
    temp=zeros(xdim,ydim,slnum,size(TE,1)); % zeros 128,128,36,5
        for i=1:size(TE,1) %1:5
        temp(:,:,:,i)=spm_read_vols(dat_echo(i)); % load image info from each TE nifti in 5x1 struct array 'echo_dat'
        end;
    dat_hdr=dat_echo; % rename header component of TE files
    dat_echo=temp; % redefine 5x1 structure array to zeros template with image info from each TE nifti 
    clear temp;

    % noise calculation
    noise=zeros(2,size(TE,1)); % zeros 2,5 (row1=mean, row2=std, column=TEs)
    ROI = spm_vol(folder_noise); % get hdr info (noise image same dimensions as TE images)
    noiseROI = spm_read_vols(ROI); % get image info

    for i=1:size(TE,1)
      dat_TE=squeeze(dat_echo(:,:,:,i)); % remove singleton dimensions from zeros template with image info from all TE niftis
      temp=noiseROI.*dat_TE; % multiply noise image info with TE nifti image info (noise ROI=1, rest zero, thus get noise image info from TE images)
      noise(1,i)=mean(nonzeros(temp(:)));
      noise(2,i)=std(nonzeros(temp(:)));
    
    end
    
    clear temp;
    cd(folder_maps)

    % initialize space for parametric maps
    
    R2s=zeros(xdim, ydim, slnum);
    N=zeros(xdim, ydim, slnum);
    N_nconv=zeros(xdim, ydim, slnum);
    S0=zeros(xdim, ydim, slnum);
    errS0=zeros(xdim, ydim, slnum);
    errR2=zeros(xdim, ydim, slnum);
    rsquare=zeros(xdim, ydim, slnum);
    rmse=zeros(xdim, ydim, slnum);

    % calculate maps
    for slice=1:slnum
       for i=1:xdim
         for j=1:ydim
              t=1; % first TE value    (Y: define good signal from each TE for each voxel)                                                            
             while t<=size(TE,1) && dat_echo(i,j,slice,t)>(5*noise(t)) % as long as TE is not the last TE and image value in pixel is 3 times noise mean at given for when t<=5 (noise(t) refernce matrix as 1 row with 10 elements, TE, then can use this voxel
                 Y(t)=dat_echo(i,j,slice,t); %a 1x? array of signal values
                 t=t+1; % repeat for voxels of next TE value until reach last TE value
             end;


                if exist('Y','var') && size(Y,2)>=10
                  X=TE(1:size(Y,2))'; %X defiine TE numbers
                  N(i,j,slice)=size(Y,2); % how many TE signal values pass noise test thus have to fit the voxel
                   
                 
              %---- create non linear fit
                    fo_ = fitoptions('method','NonlinearLeastSquares','Robust','Off','Upper',[Inf   0]);
                    st_ = [500 -20 ]; % start point for S0 and R2
                    set(fo_,'Startpoint',st_);
                    ft_ = fittype('exp1');
                    [R2starfit,goodness,O] = fit(X(:)/1000,Y(:),ft_,fo_);

                    if O.exitflag        % only consider the fits that did converge            
                    % output parameters
                        R2s(i,j,slice)=-R2starfit.b;  % use when not filter by descending order  
                        S0(i,j,slice)=R2starfit.a;
                        errors=confint(R2starfit);
                        errR2(i,j,slice)=abs(errors(1,2)-R2starfit.b)/2;
                        errS0(i,j,slice)=abs(errors(1,1)-R2starfit.a)/2;
                        rsquare(i,j,slice)=goodness.rsquare;
                        rmse(i,j,slice)=goodness.rmse;
                    else
                        N_nconv(i,j,slice)=1;
                    end
                    clear Y; % use when not filter by descending order
                end
         end
       end
    end
   

    % save maps as nifti 
    R2s_header=hdr_temp;
    R2s_header.fname=('by_voxel_R2s.nii');
    R2s_header.dt=[16 0];
    R2s_header.descrip=strcat(R2s_header.descrip,' R2s map');
    R2s_header.private.descrip=strcat(R2s_header.descrip,' R2s map');
    R2s_header.private.dat=file_array('by_voxel_R2s.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(R2s_header,R2s)
    
    S0_header=hdr_temp;
    S0_header.fname=('by_voxel_S0_star.nii');
    S0_header.dt=[16 0];
    S0_header.descrip=strcat(S0_header.descrip,' S0_star map');
    S0_header.private.descrip=strcat(S0_header.descrip,' S0_star map');
    S0_header.private.dat=file_array('by_voxel_S0_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(S0_header,S0)

    N_header=hdr_temp;
    N_header.fname=('by_voxel_N_star.nii');
    N_header.dt=[2 0];
    N_header.descrip=strcat(N_header.descrip,' number of used data points');
    N_header.private.descrip=strcat(N_header.descrip,' number of used data points');
    N_header.private.dat=file_array('by_voxel_N_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(N_header,N)
    
    N_nconv_header=hdr_temp;
    N_nconv_header.fname=('by_voxel_N_nconv_star.nii');
    N_nconv_header.dt=[2 0];
    N_nconv_header.descrip=strcat(N_nconv_header.descrip,' number of used data points');
    N_nconv_header.private.descrip=strcat(N_nconv_header.descrip,' number of used data points');
    N_nconv_header.private.dat=file_array('by_voxel_N_nconv_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(N_nconv_header,N_nconv)
    
    
    errS0_header=hdr_temp;
    errS0_header.fname=('by_voxel_errS0_star.nii');
    errS0_header.dt=[16 0];
    errS0_header.descrip=strcat(errS0_header.descrip,' errS0_star map');
    errS0_header.private.descrip=strcat(errS0_header.descrip,' errS0_star map');
    errS0_header.private.dat=file_array('by_voxel_errS0_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(errS0_header,errS0)

    errR2_header=hdr_temp;
    errR2_header.fname=('by_voxel_errR2_star.nii');
    errR2_header.dt=[16 0];
    errR2_header.descrip=strcat(errR2_header.descrip,' errR2_star map');
    errR2_header.private.descrip=strcat(errR2_header.descrip,' errR2_star map');
    errR2_header.private.dat=file_array('by_voxel_errR2_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(errR2_header,errR2)

    rsquare_header=hdr_temp;
    rsquare_header.fname=('by_voxel_rsquare_star.nii');
    rsquare_header.dt=[16 0];
    rsquare_header.descrip=strcat(rsquare_header.descrip,' rsquare_star map');
    rsquare_header.private.descrip=strcat(rsquare_header.descrip,' rsquare_star map');
    rsquare_header.private.dat=file_array('by_voxel_rsquare_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(rsquare_header,rsquare)
    
    rmse_header=hdr_temp;
    rmse_header.fname=('by_voxel_rmse_star.nii');
    rmse_header.dt=[16 0];
    rmse_header.descrip=strcat(rmse_header.descrip,' rmse_star map');
    rmse_header.private.descrip=strcat(rmse_header.descrip,' rmse_star map');
    rmse_header.private.dat=file_array('by_voxel_rmse_star.nii',[xdim ydim slnum],'FLOAT32-LE',0,1,0);
    spm_write_vol(rmse_header,rmse)
end
toc
