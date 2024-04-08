%%Script for calculating the RNA expresion levels and nascent RNA site
%%distances.

%P. Aurelio, 2019. 
%Matlab R2016a.

clear all
close all
warning off
clc

%% Some inputs:
%%----------------------------------------------------------
%Define starting directory and data types----------------------------------------
%--------------------------------------------------------------
startlocn = {"\\pasteur\SysBC-Home\schuerca\Desktop\ANALYSIS\Data\RPL3\SC\MAX"};
Condition = 'RPL3_SC';


%----------------------------------------------------------------
%Binary inputs------------------------------------------------------
%-----------------------------------------------------------------------
% Save_Results = 0; %Binary: 1 for save the results.
Save_Figures = 0; %Binary: 1 for save ALL the figures automatically in the specified folder.
Save_Results = 1; %Binary: 1 for save calculated variabled in a ,mat file.

Plot_Images = 0; %Binary: 1 for plot all the images.

subROI = 0; %1 for Analyzing subregions (This requires to input the corresponding subROI from ImageJ). If 0, this analysis won't be performed.


%----------------------------------------------------------------------------
%Some Numerical Inputs-------------------------------------------------------------
%-------------------------------------------------------------------------
Pixel_size = 65;  %Effective Pixel size in nm

PSF_roi = 5; %Roi size in pixels in order to extract the intensity from the PSFs (PSF_roi x PSF_roi). It must be an odd number.
PSF_bckg_roi = 7; %Roi size in pixels in order to take the background around the PSFs

polyScale_Init = 1; %Initial scaling factor for the polygons of around the condensates. This value is used to start the iteration.
Max_polyScale = 3; %Maximum scaling factor for the polygons of around the condensates.

Bckg_Int_thr = 1.15; %This value multiplied by the average of the minimum intensity values of each PSF of single RNAs will be used to segment the Condensate regions.

Condensates_Intensity_thr = 600; %Intensity threshold for identifying the P-Bodies and perform the colocalization step.

RNA_Condensates_dist_thr = 8; %Distance in pixels in order to recognize that an RNA localization corresponds to a Condensate. Distance between the center of mass of the condensate and the center of the RNA localization.


edges=0:0.1:12; %Histogram of Normalized Intensities edges. (Intensity/MeanIntensity)


%-----------------------------------------------------------------------------
%Initialize Variables----------------------------------------
%%----------------------------------------------------------
%--------------------------------------------------------------
PSF_length = (PSF_roi-1)/2;
PSF_bckg_length = (PSF_bckg_roi-1)/2;



%% Load the Data:


    files_RNA_images = Select1DataGroup('RNA Images','.tif',startlocn{1});
    nimages = size(files_RNA_images.data,1);

    %files_PBodies_images = Select1DataGroup('P-Bodies Images','.tif',startlocn{1});

    
% if subROI == 1;
% files_ROI{1} = Select1DataGroup('subROI (Fiji)','*.zip;*.roi',startlocn{1});
% end


%Load the files with the list of localizations and the Images
%P-Bodies images and RNA images must have the same number of frames (i.e, the same number of Z stacks)
for c=1:nimages; 
    
      if subROI == 1;
      [sROI{c}] = ReadImageJROI(strcat(files_ROI{1}.data{c,2},files_ROI{1}.data{c,1}));
      end
    

  
    %Localizations{c} = readtable(strcat(files_localizations.data{c,2},files_localizations.data{c,1})); 

    Localizations{c} = readtable(strcat(files_RNA_images.data{c,2},files_RNA_images.data{c,1}(1:end-4),'.csv'));
    Localizations_Array{c} = table2array(Localizations{c}(:,:)); 
    Localizations_Array{c}(:,2:4) = Localizations_Array{c}(:,2:4)/Pixel_size; %Transform units to pixels.
    Number_of_Localizations(c) = size(Localizations_Array{c},1); %Compute the total number of localizations for each image

%     Localizations_Array{c}(:,[2,3])= Localizations_Array{c}(:,[3,2]); %Swap X and Y columns (from ImageJ to Matlab).
    
    Im_info{c} = imfinfo(strcat(files_RNA_images.data{c,2},files_RNA_images.data{c,1}));
    numberOfFrames = length(Im_info{c});
    
            for k = 1:numberOfFrames;
                currentImage = imread(strcat(files_RNA_images.data{c,2},files_RNA_images.data{c,1}), k, 'Info', Im_info{c});
                Image_RNA{c}{k} = currentImage;
                %currentImage = imread(strcat(files_PBodies_images.data{c,2},files_PBodies_images.data{c,1}), k, 'Info', Im_info{c});
                currentImage = imread(strcat(files_RNA_images.data{c,2},files_RNA_images.data{c,1}(1:end-15),'P-Bodies_MaxProj.tif'), k, 'Info', Im_info{c});
                Image_PBodies{c}{k} = currentImage;
            end

   if Plot_Images == 1;
       figure()
       imagesc(Image_RNA{c}{1});
       title(strcat(files_RNA_images.data{c},' -- RNA Channel'),'Interpreter', 'none')
       axis square
       
       figure()
       imagesc(Image_PBodies{c}{1});
       title(strcat(files_RNA_images.data{c},' -- P-Bodies Channel'),'Interpreter', 'none')
       axis square
   end
       
end

%% Create Color variable

Color =lines(nimages);

%% Calculate some preliminary numbers

Loc2Bckg_numPX_Ratio = (PSF_roi*PSF_roi)/(PSF_bckg_roi*PSF_bckg_roi);


%% Load the subROI and extract the localizations Inside and Outside them
% if subROI == 1;
%     
%     for c=1:nimages;   
%         
%       if size(sROI{c},2) == 1;      
%      
%     [Locs_INSIDE{c}, Locs_OUTSIDE{c}] = spotsINROI(Localizations_Array{c}(:,2:3),sROI{c}); %Do the first filtering: discard nuclear localizations
%     Localizations_inROI{c} = Localizations_Array{c}(Locs_INSIDE{c}{1},:);                 
%     Localizations_outROI{c} = Localizations_Array{c}(Locs_OUTSIDE{c},:);                                                          
%                 
%     figure()
%     imagesc(Image{c}{Channel_RNA});   
%     hold on
%     scatter(sROI{c}.mnCoordinates(:,1),sROI{c}.mnCoordinates(:,2),'+');
%     hold on
%     scatter(Localizations_inROI{c}(:,2),Localizations_inROI{c}(:,3),'r');
%     hold on
%     scatter(Localizations_outROI{c}(:,2),Localizations_outROI{c}(:,3),'g');
%     title(strcat(files_images.data{c},' -- RNA Channel, ROI and segmented localizations (Green acepted; Red rejected)'));
%     axis square
%       
%       else if size(sROI{c},2) == 2;  
%               
%     [Locs_INSIDE_temp{c}, Locs_OUTSIDE_temp{c}] = spotsINROI(Localizations_Array{c}(:,2:3),sROI{c}{1}); %Do the first filtering: discard nuclear localizations
%     Localizations_inROI_temp{c} = Localizations_Array{c}(Locs_INSIDE_temp{c}{1},:);                 
%     Localizations_outROI_temp{c} = Localizations_Array{c}(Locs_OUTSIDE_temp{c},:);                                                          
%     
%     [Locs_INSIDE{c}, Locs_OUTSIDE{c}] = spotsINROI(Localizations_outROI_temp{c}(:,2:3),sROI{c}{2}); %Do the second filtering: discard localizations outside the cells
%     Localizations_inROI{c} = Localizations_outROI_temp{c}(Locs_OUTSIDE{c},:);                 
%     Localizations_outROI{c} = Localizations_outROI_temp{c}(Locs_INSIDE{c}{1},:);  
%     
%     figure()
%     imagesc(Image{c}{Channel_RNA});   
%     hold on
%     scatter(sROI{c}{1}.mnCoordinates(:,1),sROI{c}{1}.mnCoordinates(:,2),'y','+');
%     hold on
%     scatter(sROI{c}{2}.mnCoordinates(:,1),sROI{c}{2}.mnCoordinates(:,2),'y','+');
%     hold on
%     scatter(Localizations_inROI{c}(:,2),Localizations_inROI{c}(:,3),'r');
%     hold on
%     scatter(Localizations_outROI{c}(:,2),Localizations_outROI{c}(:,3),'g');
%     hold on
%     scatter(Localizations_inROI_temp{c}(:,2),Localizations_inROI_temp{c}(:,3),'r');
%     title(strcat(files_images.data{c},' -- RNA Channel, ROI and segmented localizations (Green acepted; Red rejected)'))
%     axis square         
%               
%                      
%           end
%           
%       end
%     end
%     
% else
%         for c=1:nimages;   
%    Localizations_outROI{c} = Localizations_Array{c};
%         end
%         
% end


%% Reorganize the localizations in z

  for c=1:nimages;   
      for z=1:size(Image_PBodies{c},2);
        
        idx = find((Localizations_Array{c}(:,1)==z));
        Localizations_outROI{c}{z} = Localizations_Array{c}(idx,:);
        Number_of_Localizations_perStack{c}(z) = size(Localizations_outROI{c}{z},1);

      end
  end



%% Identify the regions ocuppied by Condensates

  for c=1:nimages;   
      
      for z=1:size(Image_PBodies{c},2);
          
          Mask_idx = Image_PBodies{c}{z}>Condensates_Intensity_thr;
          Mask_idx_double = (Image_PBodies{c}{z}>Condensates_Intensity_thr).*1;
          Mask{c}{z} = double(Image_PBodies{c}{z}).*Mask_idx_double;
         
          CC{c}{z} = bwconncomp(Mask_idx);
          S{c}{z} = struct2cell(regionprops(CC{c}{z},'Centroid'));
          Condensate_Boundaries{c}{z} = bwboundaries(Mask_idx');
          Condensate_Boundaries{c}{z}(cellfun('length',Condensate_Boundaries{c}{z})<3) = []; 

          for iiii=1:size(Condensate_Boundaries{c}{z},1);
              if ~isempty(Condensate_Boundaries{c}{z});
                  Condensates_Positions{c}{z}(iiii,1) = round(mean(Condensate_Boundaries{c}{z}{iiii}(:,1)));
                  Condensates_Positions{c}{z}(iiii,2) = round(mean(Condensate_Boundaries{c}{z}{iiii}(:,2)));
              else
                Condensates_Positions{c}{z} = [];
              end
          end

          if exist('Condensates_Positions');
          else
          Condensates_Positions{c}{z} = [];
          end
            
      end

          if Plot_Images == 1;
          figure()
          imagesc(Mask{c}{1});
          title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z),' -- Mask of the Condensates'),'Interpreter', 'none');
          set(gca,'FontSize',18);
          axis square    
          end
  end

 %% If subROIs, remove Condensates that are outside the specified region
 
 if subROI == 1;
    
    for c=1:nimages;   
        
      if size(sROI{c},2) == 1;      
     
     [Cond_INSIDE{c}, Cond_OUTSIDE{c}] = spotsINROI(Condensates_Positions{c},sROI{c}); %Do the first filtering: discard nuclear localizations
    Cond_inROI{c} = Condensates_Positions{c}(Cond_INSIDE{c}{1},:);                 
    Cond_outROI{c} = Condensates_Positions{c}(Cond_OUTSIDE{c},:);                                                          
      
      else if size(sROI{c},2) == 2;  
              
    [Cond_INSIDE_temp{c}, Cond_OUTSIDE_temp{c}] = spotsINROI(Condensates_Positions{c},sROI{c}{1}); %Do the first filtering: discard nuclear localizations
    Cond_inROI_temp{c} = Condensates_Positions{c}(Cond_INSIDE_temp{c}{1},:);                 
    Cond_outROI_temp{c} = Condensates_Positions{c}(Cond_OUTSIDE_temp{c},:);                                                          
    
    [Cond_INSIDE{c}, Cond_OUTSIDE{c}] = spotsINROI(Cond_outROI_temp{c},sROI{c}{2}); %Do the second filtering: discard localizations outside the cells
    Cond_inROI{c} = Cond_outROI_temp{c}(Cond_OUTSIDE{c},:);                 
    Cond_outROI{c} = Cond_outROI_temp{c}(Cond_INSIDE{c}{1},:);  
           
                     
          end
          
      end
    end
    
 else

     Cond_outROI = Condensates_Positions;

end

  
  
%% Perform Colocalization and remove the localizations corresponding to condensates from the list of RNA localizations

  for c=1:nimages;   
      for z=1:size(Image_PBodies{c},2);

       Localizations_outROI_all{c}{z} = Localizations_outROI{c}{z};
    

     if ~isempty(Localizations_outROI{c}{z}(:,2:3)) && ~isempty(Condensate_Boundaries{c}{z});
    
    [Overlapping_locs_boundaries{c}{z}, othr,Positive_Cond_idx{c}{z}] = spotsINROI_V3(Localizations_outROI{c}{z}(:,2:3), Condensate_Boundaries{c}{z});   %Identify Localizations inside Condensate Boundaries

    Positive_Condensate_list{c}{z} = unique(nonzeros(Positive_Cond_idx{c}{z}));
    
    Distance_RNA_Cond{c}{z} = pdist2(Cond_outROI{c}{z},Localizations_outROI{c}{z}(:,2:3))'; %Identify Localizations that are closer than the specified distance to the center of the condensates
     [Overlapping_locs_idx_dist{c}{z}, othr] = find(Distance_RNA_Cond{c}{z} < RNA_Condensates_dist_thr);
 
    Overlapping_locs_idx{c}{z} = [vertcat(Overlapping_locs_boundaries{c}{z}{:}); Overlapping_locs_idx_dist{c}{z}];
    Overlapping_locs_idx{c}{z} = unique(Overlapping_locs_idx{c}{z});

  
    
  if Plot_Images == 1;
  figure()
  imagesc(Image_RNA{c}{z});
  hold on
  scatter(Localizations_outROI{c}{z}(:,2),Localizations_outROI{c}{z}(:,3),'g');
  hold on
  scatter(Cond_outROI{c}{z}(:,1),Cond_outROI{c}{z}(:,2),'y','filled');
  hold on
  scatter(Localizations_outROI{c}{z}(Overlapping_locs_idx{c}{z},2),Localizations_outROI{c}{z}(Overlapping_locs_idx{c}{z},3),'r');
  hold on
  for iii=1:size(Condensate_Boundaries{c}{z},1);
  scatter(Condensate_Boundaries{c}{z}{iii}(:,1),Condensate_Boundaries{c}{z}{iii}(:,2),3,'y','+');
  hold on   
  end
  axis square
    
  aa = [1:size(Cond_outROI{c}{z},1)]'; 
  bb = num2str(aa); 
  cc = cellstr(bb);
  labelpoints(Cond_outROI{c}{z}(:,1),Cond_outROI{c}{z}(:,2),cc,'S',0.3)


%   dx = 0.1; dy = 0.1;
%   text(Condensates_Positions{c}(:,1)+dx, Condensates_Positions{c}(:,2)+dy, c);
  title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z),' -- RNA Localizations (green); Condensate Boundaries and Center Positions (yellow); Overlapping RNA Localizations (red)'),'Interpreter', 'none')
  set(gca,'FontSize',18);
  end



  Overlapping_Localizations{c}{z} = Localizations_outROI{c}{z}(Overlapping_locs_idx{c}{z},:);

  Localizations_outROI{c}{z}(Overlapping_locs_idx{c}{z},:) = [];

  Total_number_localizations{c}(z) = size(Localizations_outROI_all{c}{z},1);
  Total_number_Overlapping_localizations{c}(z) = size(Overlapping_Localizations{c}{z},1);


  Total_number_condensates{c}(z) = size(Cond_outROI{c}{z},1);
  Total_number_Positive_condensates{c}(z) = size(Positive_Condensate_list{c}{z},1);

     else
     
   Total_number_localizations{c}(z) = nan;
   Total_number_Overlapping_localizations{c}(z) = nan;

  Total_number_condensates{c}(z) = nan;
  Total_number_Positive_condensates{c}(z) = nan;

  Positive_Condensate_list{c}{z} = [];

     end
       
      end
    
      RNA_Colocalization_Percentage(c) = nansum(Total_number_Overlapping_localizations{c})*100/nansum(Total_number_localizations{c});

      Condensate_Positive_Percentage(c) = nansum(Total_number_Positive_condensates{c})*100/nansum(Total_number_condensates{c});

  end



 %% Print Colozalization Results


  for c=1:nimages;   

  fprintf('\n');    
    fprintf('---------------------------------');
    fprintf('\n'); 
    fprintf(strcat('Image_',num2str(c)));   
    fprintf('\n');    
    fprintf(files_RNA_images.data{c});  
    fprintf('\n');   
    fprintf('Total percentage of RNAS inside condensates:\n');
    fprintf('------');
    fprintf('\n'); 
    fprintf(strcat(num2str(RNA_Colocalization_Percentage(c)),' %%'));     
    fprintf('\n');
    fprintf('\n'); 
    fprintf('Total percentage of positive Condensates:\n');   
    fprintf('------');
    fprintf('\n');
    fprintf(strcat(num2str(Condensate_Positive_Percentage(c)),' %%'));
    fprintf('\n');
    fprintf('---------------------------------');
    fprintf('\n');      
    fprintf('\n');    
end


%% Obtain the Intensity values of each localization

for c=1:nimages;

      for z=1:size(Image_PBodies{c},2);

        
        
        PSF_x{c}{z} = round(Localizations_outROI{c}{z}(:,3));
        PSF_y{c}{z} = round(Localizations_outROI{c}{z}(:,2));
         
                for i=1:length(PSF_x{c}{z});
                    
                    if (PSF_x{c}{z}(i) - PSF_bckg_length) > 0 & (PSF_x{c}{z}(i) + PSF_bckg_length) <= size(Image_RNA{c}{z},1) & (PSF_y{c}{z}(i) - PSF_bckg_length) > 0 & (PSF_y{c}{z}(i) + PSF_bckg_length) <= size(Image_RNA{c}{z},2);
                 PSF_Int{c}{z}(i) = sum(sum(Image_RNA{c}{z}(PSF_x{c}{z}(i) - PSF_length:PSF_x{c}{z}(i) + PSF_length, PSF_y{c}{z}(i) - PSF_length: PSF_y{c}{z}(i) + PSF_length)));
                 PSF_Bckg{c}{z}(i) = sum(sum(Image_RNA{c}{z}(PSF_x{c}{z}(i) - PSF_bckg_length:PSF_x{c}{z}(i) + PSF_bckg_length, PSF_y{c}{z}(i) - PSF_bckg_length: PSF_y{c}{z}(i) + PSF_bckg_length)));
               
                 I1 = PSF_Int{c}{z}(i);
                 I2 = PSF_Bckg{c}{z}(i);
                 A1 = PSF_roi*PSF_roi;
                 A2 = PSF_bckg_roi*PSF_bckg_roi;
                 
                 bckg_level{c}{z}(i) = min(min(Image_RNA{c}{z}(PSF_x{c}{z}(i) - PSF_bckg_length:PSF_x{c}{z}(i) + PSF_bckg_length, PSF_y{c}{z}(i) - PSF_bckg_length: PSF_y{c}{z}(i) + PSF_bckg_length)));
                 
                Localization_Normalized_Intensity{c}{z}(i) = I1 - ((I2-I1)*A1/(A2-A1));
                
                    else
                
                Localization_Normalized_Intensity{c}{z}(i) = 0;

                    end               
                    
                end
                             
end

end
      
%% Fit the CDF to the Intensities distribution

                        for c=1:nimages; 

                                  for z=1:size(Image_PBodies{c},2);
                       

                                Localization_Normalized_Intensity_all{c}{z} = nonzeros(Localization_Normalized_Intensity{c}{z}');
                                Localization_Normalized_Intensity_all{c}{z} = Localization_Normalized_Intensity_all{c}{z}(Localization_Normalized_Intensity_all{c}{z} > 0);
                                Localization_Normalized_Intensity_all{c}{z} = rmoutliers(Localization_Normalized_Intensity_all{c}{z});
                                

                              f2_all{c}{z} = fitdist(Localization_Normalized_Intensity_all{c}{z},'Normal');
                              f4_all{c}{z} = fitdist(Localization_Normalized_Intensity_all{c}{z},'Exp');
                              Mean_Brightness_Single_RNA{c}{z} = f2_all{c}{z}.mu;
%                              Std_Brightness_Single_RNA{c} = std(Localization_Normalized_Intensity_all{c});
                              ci2_all{c}{z} = paramci(f2_all{c}{z}); %Get 95% Confidence Interval
                              [f1_all,x1_all] = ecdf(Localization_Normalized_Intensity_all{c}{z});
                                
                                
                               if Plot_Images == 1;
                                figure()
                                scatter(x1_all,f1_all,'filled');
                                hold on
                                f3_all{c}{z} = cdf(f2_all{c}{z},x1_all);
                                plot(x1_all,f3_all{c}{z});
                                title(strcat('Image_',num2str(c),' - Z= ',num2str(z)),'Interpreter','none');
                                xlabel('Intensity');
                                ylabel('CDF');
                                title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z)),' -- Cumulative Density Function of RNA Intensities (Normal Dist Fit)','Interpreter','none')
                                set(gca,'FontSize',18);
              
                                
%                                figure()
%                                 scatter(x1_all,f1_all,'MarkerFaceColor',Color(c,:));
%                                 hold on
%                                 f5_all{c} = cdf(f4_all{c},x1_all);
%                                 plot(x1_all,f5_all{c});
%                                 title(strcat('Image_',num2str(c)),'Interpreter','none');
%                                 xlabel('Intensity');
%                                 ylabel('CDF');
%                                 hold on
%                                 title('Cumulative Density Function of Localization Intensities (Exponential Dist Fit)');   
                                
                                figure()
                                hist(Localization_Normalized_Intensity_all{c}{z},25);
                                title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z)),' -- Histogram of RNA Intensities','Interpreter','none');
                                 xlabel('Intensities');
                                 ylabel('Counts');
                                 set(gca,'FontSize',18);

                               end
                                
    fprintf('\n');    
    fprintf('---------------------------------');
    fprintf('\n'); 
    fprintf(strcat('Image_',num2str(c)));   
    fprintf('\n');    
    fprintf(files_RNA_images.data{c});  
    fprintf('\n'); 
    fprintf(strcat('Z plane = ',num2str(z)));   
    fprintf('\n');    
    fprintf('Single RNA mean intensity (ADU Counts):\n');
    fprintf(num2str(round(Mean_Brightness_Single_RNA{c}{z})));
    fprintf('\n');
    fprintf('Single RNA mean intensity range - 95%% Confidence Interval (ADU Counts):\n');
    fprintf(strcat('[',num2str(round(ci2_all{c}{z}(1,1))),'-',num2str(round(ci2_all{c}{z}(2,1))),']'));
    fprintf('\n');
    fprintf('---------------------------------');
    fprintf('\n');
                         
                         
  end   

  end

%% Plot the distribution of intensities
                       

%                         for c=1:nimages; 
%                       
%                                   for z=1:size(Image_PBodies{c},2);
%                                         
%                                       if Plot_Images == 1;
%                                 figure()
%                                 scatter(1:1:length(Localization_Normalized_Intensity_all{c}{z}),Localization_Normalized_Intensity_all{c}{z});
%                                 title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z),' -- Normalized intensities of the localizations'),'Interpreter','none');
%                                 xlabel('Localization');
%                                 ylabel('Intensity');
%                                 set(gca,'FontSize',18);
%                                       end
%                       
%                                 end
% 
%                         end
                        
  %% Calculate the intensity of each Condensate and the number of RNAs on each Condensate

  polyScale = polyScale_Init;



    for c=1:nimages;

            for z=1:size(Image_PBodies{c},2);
                    
                if ~isempty(Condensate_Boundaries{c}{z});

                PSF_Cond_x{c}{z} = round(Cond_outROI{c}{z}(:,2));
                PSF_Cond_y{c}{z} = round(Cond_outROI{c}{z}(:,1));
        
                for i=1:length(PSF_Cond_x{c}{z});
                                                
                            polyin{c}{z}{i} = polyshape(Condensate_Boundaries{c}{z}{i});
                            poly2{c}{z}{i} = scale(polyin{c}{z}{i},polyScale,[Cond_outROI{c}{z}(i,1) Cond_outROI{c}{z}(i,2)]);
                            polyMASK{c}{z}{i} = poly2mask(rmmissing(poly2{c}{z}{i}.Vertices(:,1)),rmmissing(poly2{c}{z}{i}.Vertices(:,2)),size(Image_RNA{c}{z},1),size(Image_RNA{c}{z},2));
                            polyPIXELS{c}{z}{i} = find(polyMASK{c}{z}{i}>0); 
                             
                            if length(polyPIXELS{c}{z}{i})>9;
            
                             PSF_Cond_Int_Region{c}{z}{i} = Image_RNA{c}{z}(polyPIXELS{c}{z}{i});
                             idx_Int{c}{z}{i} = find(PSF_Cond_Int_Region{c}{z}{i} > Bckg_Int_thr*mean(bckg_level{c}{z}));
                             
                             ratio_diff =  (length(idx_Int{c}{z}{i})/length(PSF_Cond_Int_Region{c}{z}{i}))-Loc2Bckg_numPX_Ratio;
                             
                            while (ratio_diff>0.05 && polyScale<Max_polyScale);      %Iterate in order to achieve a ratio between the size of the background and the size of the Condensate region that is equivalent to the one specified for the RNA localizations
                           
                             poly2{c}{z}{i} = scale(polyin{c}{z}{i},polyScale,[Cond_outROI{c}{z}(i,1) Cond_outROI{c}{z}(i,2)]);
                             polyMASK{c}{z}{i} = poly2mask(rmmissing(poly2{c}{z}{i}.Vertices(:,1)),rmmissing(poly2{c}{z}{i}.Vertices(:,2)),size(Image_RNA{c}{z},1),size(Image_RNA{c}{z},2));
                             polyPIXELS{c}{z}{i} = find(polyMASK{c}{z}{i}>0);
            
                             PSF_Cond_Int_Region{c}{z}{i} = Image_RNA{c}{z}(polyPIXELS{c}{z}{i});
                             idx_Int{c}{z}{i} = find(PSF_Cond_Int_Region{c}{z}{i} > Bckg_Int_thr*mean(bckg_level{c}{z}));            
                           
                             ratio_diff =  abs((length(idx_Int{c}{z}{i})/length(PSF_Cond_Int_Region{c}{z}{i}))-Loc2Bckg_numPX_Ratio);      
            
                             polyScale = polyScale + 0.05; 
                            end
                            
                             PSF_Cond_Int{c}{z}(i) = sum(sum(PSF_Cond_Int_Region{c}{z}{i}(idx_Int{c}{z}{i})));
                             PSF_Cond_Bckg{c}{z}(i) = sum(sum(PSF_Cond_Int_Region{c}{z}{i}));
                           
                             I1_Cond = PSF_Cond_Int{c}{z}(i);
                             I2_Cond = PSF_Cond_Bckg{c}{z}(i);
                             A1_Cond = length(idx_Int{c}{z}{i});
                             A2_Cond = length(PSF_Cond_Int_Region{c}{z}{i});
                             
                            Localization_Normalized_Intensity_Cond{c}{z}(i) = I1_Cond - ((I2_Cond-I1_Cond)*A1_Cond/(A2_Cond-A1_Cond));
                        
                            Number_of_RNAS_per_Condensate{c}{z}(i) = Localization_Normalized_Intensity_Cond{c}{z}(i)/Mean_Brightness_Single_RNA{c}{z};
                            Std_Sup_Number_of_RNAS_per_Condensate{c}{z}(i) = Localization_Normalized_Intensity_Cond{c}{z}(i)/ci2_all{c}{z}(1,1);
                            Std_Inf_Number_of_RNAS_per_Condensate{c}{z}(i) = Localization_Normalized_Intensity_Cond{c}{z}(i)/ci2_all{c}{z}(2,1);
                            
                            polyScale = polyScale_Init; %Re-initialize variable


                            else

                               Localization_Normalized_Intensity_Cond{c}{z}(i) = 0;

                               Localization_Normalized_Intensity_Cond{c}{z}(i) = 0;
                               Number_of_RNAS_per_Condensate{c}{z}(i) = 0;
                               Std_Sup_Number_of_RNAS_per_Condensate{c}{z}(i) = 0;
                               Std_Inf_Number_of_RNAS_per_Condensate{c}{z}(i) = 0;
                                            
                            end                     
                            


                end

                            
               Percentage_of_RNAS_inside_Condensates{c}{z} = sum(Number_of_RNAS_per_Condensate{c}{z},'omitnan')/(sum(Number_of_RNAS_per_Condensate{c}{z},'omitnan') + size(Localizations_outROI{c}{z},1));
               Std_Sup_Percentage_of_RNAS_inside_Condensates{c}{z} = sum(Std_Sup_Number_of_RNAS_per_Condensate{c}{z},'omitnan')/(sum(Std_Sup_Number_of_RNAS_per_Condensate{c}{z},'omitnan') + size(Localizations_outROI{c}{z},1));
               Std_Inf_Percentage_of_RNAS_inside_Condensates{c}{z} = sum(Std_Inf_Number_of_RNAS_per_Condensate{c}{z},'omitnan')/(sum(Std_Inf_Number_of_RNAS_per_Condensate{c}{z},'omitnan') + size(Localizations_outROI{c}{z},1));
             



                else

                   Percentage_of_RNAS_inside_Condensates{c}{z} = [];
                   Std_Sup_Percentage_of_RNAS_inside_Condensates{c}{z} = [];
                   Std_Inf_Percentage_of_RNAS_inside_Condensates{c}{z} = [];

                   Localization_Normalized_Intensity_Cond{c}{z} = [];

                   Number_of_RNAS_per_Condensate{c}{z} = [];
                   Std_Sup_Number_of_RNAS_per_Condensate{c}{z} = [];
                   Std_Inf_Number_of_RNAS_per_Condensate{c}{z} = [];
                            

                end



               Localization_Normalized_Intensity_Cond_Positive{c}{z} = nonzeros(Localization_Normalized_Intensity_Cond{c}{z});





               
    fprintf('\n');    
    fprintf('---------------------------------');
    fprintf('\n'); 
    fprintf(strcat('Image_',num2str(c)));   
    fprintf('\n');    
    fprintf(files_RNA_images.data{c});  
    fprintf('\n'); 
    fprintf(strcat('Z plane = ',num2str(z)));   
    fprintf('\n');    
    fprintf('Number of RNAS inside condensates:\n');
    fprintf('------');
    fprintf('\n'); 
    
    for ii = 1:size(Number_of_RNAS_per_Condensate{c}{z},2);

        if ~isempty(find(Positive_Condensate_list{c}{z}==ii));
            fprintf(strcat('Condensate_',num2str(ii),' -- Mean: ',num2str(Number_of_RNAS_per_Condensate{c}{z}(ii)),...
            ' -- Range (95%%CI): [',num2str(Std_Inf_Number_of_RNAS_per_Condensate{c}{z}(ii)),'-',num2str(Std_Sup_Number_of_RNAS_per_Condensate{c}{z}(ii)),']',...
            ' -- Region size (px): ',num2str(length(PSF_Cond_Int_Region{c}{z}{ii})),' -- Condensate size (px): ',num2str(length(PSF_Cond_Int_Region{c}{z}{ii}(idx_Int{c}{z}{ii}))),...
            ' -- Density (RNAs/px): ',num2str(Number_of_RNAS_per_Condensate{c}{z}(ii)/length(PSF_Cond_Int_Region{c}{z}{ii}(idx_Int{c}{z}{ii})))));     
            fprintf('\n');

        else
            fprintf(strcat('Condensate_',num2str(ii),' -- Mean: ',num2str(Number_of_RNAS_per_Condensate{c}{z}(ii)),...
            ' -- Range (95%%CI): [',num2str(Std_Inf_Number_of_RNAS_per_Condensate{c}{z}(ii)),'-',num2str(Std_Sup_Number_of_RNAS_per_Condensate{c}{z}(ii)),']'));
            fprintf('\n');
        end

    end

    fprintf('\n'); 
    fprintf('Percentage of RNAS inside condensates for each Z Plane:\n');   
    fprintf('------');
    fprintf('\n');
    fprintf(strcat('Image_',num2str(c)));   
    fprintf('\n');    
    fprintf(files_RNA_images.data{c});  
    fprintf('\n'); 
    fprintf(strcat('Z plane = ',num2str(z)));   
    fprintf('\n');    
    fprintf(strcat('Mean: ',num2str(Percentage_of_RNAS_inside_Condensates{c}{z}*100),'%% ',...
        ' -- Range: [',num2str(Std_Inf_Percentage_of_RNAS_inside_Condensates{c}{z}*100),'-',num2str(Std_Sup_Percentage_of_RNAS_inside_Condensates{c}{z}*100),']%%'));
    fprintf('\n');
    fprintf('---------------------------------');
    fprintf('\n');      
    fprintf('\n');    
                  
  
  if Plot_Images == 1;  
  figure()
  imagesc(Image_RNA{c}{z});
  hold on

  for iii=1:size(PSF_Cond_Int_Region{c}{z},2);
  
  if ~isempty(polyin{c}{z}{iii}.Vertices);    
  polyin_plot{c}{z}{iii}.Vertices(:,:) = [polyin{c}{z}{iii}.Vertices(:,:); polyin{c}{z}{iii}.Vertices(1,:)];
  else
  polyin_plot{c}{z}{iii}.Vertices = polyin{c}{z}{iii}.Vertices;    
  end
  plot(polyin_plot{c}{z}{iii}.Vertices(:,1),polyin_plot{c}{z}{iii}.Vertices(:,2),'g','LineWidth',1);
  hold on

  
  if ~isempty(poly2{c}{z}{iii}.Vertices);    
  poly2_plot{c}{z}{iii}.Vertices(:,:) = [poly2{c}{z}{iii}.Vertices(:,:); poly2{c}{z}{iii}.Vertices(1,:)];
  else
  poly2_plot{c}{z}{iii}.Vertices = poly2{c}{z}{iii}.Vertices;
  end
  plot(poly2_plot{c}{z}{iii}.Vertices(:,1),poly2_plot{c}{z}{iii}.Vertices(:,2),'y','LineWidth',1);
  hold on
  

  [row,column] = ind2sub(size(Image_RNA{c}{z}),polyPIXELS{c}{z}{iii});
  scatter(column,row,5,'+','k');
  hold on
  [row,column] = ind2sub(size(Image_RNA{c}{z}),polyPIXELS{c}{z}{iii}(idx_Int{c}{z}{iii}));
  scatter(column,row,5,'+','w');
  hold on
 

  end

  aa = [1:size(PSF_Cond_Int_Region{c}{z},2)]'; 
  bb = num2str(aa); 
  cc = cellstr(bb);
  labelpoints(Cond_outROI{c}{z}(1:size(PSF_Cond_Int_Region{c}{z},2),1),Cond_outROI{c}{z}(1:size(PSF_Cond_Int_Region{c}{z},2),2),cc,'S',0.1);
    
  axis square
  colormap default
  title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z),' -- Original (green) and Scaled (yellow) Condensate Regions, Condensate Pixels (white), Condensate Background Pixels (black)'),'Interpreter','None')
    set(gca,'FontSize',18);
    
  


    figure()
  imagesc(Image_PBodies{c}{z});
  hold on

  for iii=1:size(PSF_Cond_Int_Region{c}{z},2);
  
  if ~isempty(polyin{c}{z}{iii}.Vertices);    
  polyin_plot{c}{z}{iii}.Vertices(:,:) = [polyin{c}{z}{iii}.Vertices(:,:); polyin{c}{z}{iii}.Vertices(1,:)];
  else
  polyin_plot{c}{z}{iii}.Vertices = polyin{c}{z}{iii}.Vertices;    
  end
  plot(polyin_plot{c}{z}{iii}.Vertices(:,1),polyin_plot{c}{z}{iii}.Vertices(:,2),'g','LineWidth',1);
  hold on

  
  if ~isempty(poly2{c}{z}{iii}.Vertices);    
  poly2_plot{c}{z}{iii}.Vertices(:,:) = [poly2{c}{z}{iii}.Vertices(:,:); poly2{c}{z}{iii}.Vertices(1,:)];
  else
  poly2_plot{c}{z}{iii}.Vertices = poly2{c}{z}{iii}.Vertices;
  end
  plot(poly2_plot{c}{z}{iii}.Vertices(:,1),poly2_plot{c}{z}{iii}.Vertices(:,2),'y','LineWidth',1);
  hold on
  

  [row,column] = ind2sub(size(Image_RNA{c}{z}),polyPIXELS{c}{z}{iii});
  scatter(column,row,5,'+','k');
  hold on
  [row,column] = ind2sub(size(Image_RNA{c}{z}),polyPIXELS{c}{z}{iii}(idx_Int{c}{z}{iii}));
  scatter(column,row,5,'+','w');
  hold on
 

  end

  aa = [1:size(PSF_Cond_Int_Region{c}{z},2)]'; 
  bb = num2str(aa); 
  cc = cellstr(bb);
  labelpoints(Cond_outROI{c}{z}(1:size(PSF_Cond_Int_Region{c}{z},2),1),Cond_outROI{c}{z}(1:size(PSF_Cond_Int_Region{c}{z},2),2),cc,'S',0.1);
    
  axis square
  colormap default
  title(strcat(files_RNA_images.data{c},' - Z= ',num2str(z),' -- Original (green) and Scaled (yellow) Condensate Regions, Condensate Pixels (white), Condensate Background Pixels (black)'),'Interpreter','None')
    set(gca,'FontSize',18);  

  end


               

            
  end
    






    Percentage_of_RNAS_inside_Condensates_Mean{c} = mean([Percentage_of_RNAS_inside_Condensates{c}{:}]);
    Std_Sup_Percentage_of_RNAS_inside_Condensates_Mean{c} = mean([Std_Sup_Percentage_of_RNAS_inside_Condensates{c}{:}]);
    Std_Inf_Percentage_of_RNAS_inside_Condensates_Mean{c} = mean([Std_Inf_Percentage_of_RNAS_inside_Condensates{c}{:}]);


fprintf('\n'); 
    fprintf('Overall percentage of RNAS inside condensates for each Image:\n');   
    fprintf('------');
    fprintf('\n');
    fprintf(strcat('Image_',num2str(c)));   
    fprintf('\n');    
    fprintf(files_RNA_images.data{c});   
    fprintf('\n');    
    fprintf(strcat('Mean: ',num2str(Percentage_of_RNAS_inside_Condensates_Mean{c}*100),'%% ',...
        ' -- Range: [',num2str(Std_Inf_Percentage_of_RNAS_inside_Condensates_Mean{c}*100),'-',num2str(Std_Sup_Percentage_of_RNAS_inside_Condensates_Mean{c}*100),']%%'));
    fprintf('\n');
    fprintf('---------------------------------');
    fprintf('\n');      
    fprintf('\n');    
                  


    end






%% Calculate Normalize Intensity and Plot Histograms of Intensity of RNA in Condensates and Intensity of RNA outside
%%Condensates

                       figure()

                        for c=1:nimages; 

                                for z=1:size(Image_PBodies{c},2);

                                Localization_Normalized_Intensity_RNAs{c}{z} = Localization_Normalized_Intensity_all{c}{z}/Mean_Brightness_Single_RNA{c}{z};
                                Localization_Normalized_Intensity_Condensates{c}{z} = Localization_Normalized_Intensity_Cond{c}{z}(Positive_Condensate_list{c}{z})/Mean_Brightness_Single_RNA{c}{z}; %Only Positive Condensates!!
                                Localization_Normalized_Intensity_Condensates{c}{z} = Localization_Normalized_Intensity_Condensates{c}{z}';

                                Localization_Intensity_RNAs{c}{z} = Localization_Normalized_Intensity_all{c}{z};
                                Localization_Intensity_Condensates{c}{z} = Localization_Normalized_Intensity_Cond{c}{z}(Positive_Condensate_list{c}{z});

                                [N,X] = hist(Localization_Normalized_Intensity_RNAs{c}{z},edges);
                                N = N/sum(N); %Normalize Histogram.Sum of counts equals 1.
                                bar(X,N,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.6);
                                hold on
                                [N,X] = hist(Localization_Normalized_Intensity_Condensates{c}{z},edges);
                                N = N/sum(N); %Normalize Histogram.Sum of counts equals 1.
                                bar(X,N,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.6);
%                                 xlim([])
%                                 ylim([])

                                  end
                        end
                        
                        title(strcat(files_RNA_images.data{c},' -- Histogram of Normalized RNA Intensities'),'Interpreter','None');
                        xlabel('Normalized Intensities');
                        ylabel('Counts');
                        legend({'RNA Outside','RNA in P-Bodies'});
                        set(gca,'FontSize',18);




 %% Collapse all images and plot Histograms of Intensity of RNA in Condensates and Intensity of RNA outside


     for c=1:nimages; 

         Localization_Normalized_Intensity_RNAs_Collapse_images{c} = cat(1,Localization_Normalized_Intensity_RNAs{c}{:});
         Localization_Normalized_Intensity_Condensates_Collapse_images{c} = cat(1,Localization_Normalized_Intensity_Condensates{c}{:});

         Localization_Intensity_RNAs_Collapse_images{c} = cat(1,Localization_Intensity_RNAs{c}{:});
         Localization_Intensity_Condensates_Collapse_images{c} =  cat(1,Localization_Intensity_Condensates{c}{:});

         Percentage_of_RNAS_inside_Condensates_Collapse_images{c} = cat(1,Percentage_of_RNAS_inside_Condensates{c}{:});
         Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_images{c} = cat(1,Std_Inf_Percentage_of_RNAS_inside_Condensates{c}{:});
         Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_images{c} = cat(1,Std_Sup_Percentage_of_RNAS_inside_Condensates{c}{:});

     end


      Localization_Normalized_Intensity_RNAs_Collapse_all = cat(1,Localization_Normalized_Intensity_RNAs_Collapse_images{:});
      Localization_Intensity_RNAs_Collapse_all = cat(1,Localization_Intensity_RNAs_Collapse_images{:});

      emptyCells = cellfun(@isempty,Localization_Normalized_Intensity_Condensates_Collapse_images);
      Localization_Normalized_Intensity_Condensates_Collapse_images(emptyCells) = {0};
      Localization_Normalized_Intensity_Condensates_Collapse_all = cat(1,Localization_Normalized_Intensity_Condensates_Collapse_images{:});      

      Localization_Intensity_Condensates_Collapse_images(cellfun(@isempty,Localization_Intensity_Condensates_Collapse_images)) = {0};
      Localization_Intensity_Condensates_Collapse_all = cat(1,Localization_Intensity_Condensates_Collapse_images{:});

      %Remove NaNs
      Localization_Normalized_Intensity_RNAs_Collapse_all(isnan(Localization_Normalized_Intensity_RNAs_Collapse_all)) = [];
      Localization_Normalized_Intensity_Condensates_Collapse_all(isnan(Localization_Normalized_Intensity_Condensates_Collapse_all)) = [];
      Localization_Intensity_RNAs_Collapse_all(isnan(Localization_Intensity_RNAs_Collapse_all)) = [];
      Localization_Intensity_Condensates_Collapse_all(isnan(Localization_Intensity_Condensates_Collapse_all)) = [];
      %Remove Zeros
      Localization_Normalized_Intensity_RNAs_Collapse_all = nonzeros(Localization_Normalized_Intensity_RNAs_Collapse_all);
      Localization_Normalized_Intensity_Condensates_Collapse_all = nonzeros(Localization_Normalized_Intensity_Condensates_Collapse_all);
      Localization_Intensity_RNAs_Collapse_all = nonzeros(Localization_Normalized_Intensity_RNAs_Collapse_all);
      Localization_Intensity_Condensates_Collapse_all = nonzeros(Localization_Intensity_Condensates_Collapse_all);

      Percentage_of_RNAS_inside_Condensates_Collapse_all = cat(1,Percentage_of_RNAS_inside_Condensates_Collapse_images{:});
      Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_all = cat(1,Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_images{:});
      Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_all = cat(1,Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_images{:});


        figure()

        [N,X] = hist(Localization_Normalized_Intensity_RNAs_Collapse_all,edges);
        N = N/sum(N); %Normalize Histogram.Sum of counts equals 1.
        bar(X,N,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.6);
        hold on
        [N,X] = hist(Localization_Normalized_Intensity_Condensates_Collapse_all,edges);
        N = N/sum(N); %Normalize Histogram.Sum of counts equals 1.
        bar(X,N,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.6);
        
        title(strcat(files_RNA_images.data{c},' -- Histogram of Normalized RNA Intensities'),'Interpreter','None');
        xlabel('Normalized Intensities');
        ylabel('Counts');
        legend({'RNA Outside','RNA in P-Bodies'});
        set(gca,'FontSize',18);
%                                 xlim([])
%                                 ylim([])

    

%% Save figures and results


if Save_Results == 1
    
      for c=1:nimages;

            for z=1:size(Image_PBodies{c},2);

                for ii = 1:size(Number_of_RNAS_per_Condensate{c}{z},2);

                    if ~isempty(PSF_Cond_Int_Region{c});
                            if ~isempty(PSF_Cond_Int_Region{c}{z});
                                if ii>size(PSF_Cond_Int_Region{c}{z},2);
                                        RNAs_Density_in_Condensates{c}{z}(ii) = 0;    
                                        Condensate_Region_Size{c}{z}(ii) = 0;
                                        Condensate_Int_Region_Size{c}{z}(ii) = 0;
                                else

                                        RNAs_Density_in_Condensates{c}{z}(ii) = Number_of_RNAS_per_Condensate{c}{z}(ii)/length(PSF_Cond_Int_Region{c}{z}{ii}(idx_Int{c}{z}{ii}));    
                                        Condensate_Region_Size{c}{z}(ii) = length(PSF_Cond_Int_Region{c}{z}{ii});
                                        Condensate_Int_Region_Size{c}{z}(ii) = length(PSF_Cond_Int_Region{c}{z}{ii}(idx_Int{c}{z}{ii}));

                                end

                        else
                        RNAs_Density_in_Condensates{c}{z}(ii) = 0;    
                        Condensate_Region_Size{c}{z}(ii) = 0;
                        Condensate_Int_Region_Size{c}{z}(ii) = 0;

                        end

                    else
                        RNAs_Density_in_Condensates{c}{z}(ii) = 0;    
                        Condensate_Region_Size{c}{z}(ii) = 0;
                        Condensate_Int_Region_Size{c}{z}(ii) = 0;

                    end


                end

            end
                
        Results{c}{1,1} = 'Image name';
        Results{c}{1,2} = 'Mean intensity of single RNA [ADU]';
        Results{c}{1,3} = 'Mean intensity of single RNA 95 Confidence Intervals [ADU]';
        Results{c}{1,4} = 'Number of RNAs per Condensate';
        Results{c}{1,5} = 'Number of RNAs per Condensate Inferior 95 Confidence Interval';
        Results{c}{1,6} = 'Number of RNAs per Condensate Superior 95 Confidence Interval';
        Results{c}{1,7} = 'Total Condensate Region Size (Condensate and Background) [pixel]';
        Results{c}{1,8} = 'Condensate Region Size (Intensity) [pixel]';
        Results{c}{1,9} =  'RNAs Density in Condensates [RNAs/pixel]';
        Results{c}{1,10} = 'Percentage of RNA Inside Condensates (0 to 1)';
        Results{c}{1,11} =  'Percentage of RNA Inside Condensates Inferior 95 Confidence Interval (0 to 1)';
        Results{c}{1,12} = 'Percentage of RNA Inside Condensates Superior 95 Confidence Interval (0 to 1)';       
        Results{c}{1,13} = 'RNAs Intensities (Outside Condensates)';     
        Results{c}{1,14} = 'RNAs Intensities (Inside Condensates) (Only Positive Condensates)'; 
        Results{c}{1,15} = 'Total number of localizations per image'; 
        Results{c}{1,16} = 'Number of localizations per image per z-stack'; 
                
        if exist('files_RNA_images')
        Results{c}{2,1} = files_RNA_images.data{c};
        end
        if exist('Mean_Brightness_Single_RNA')
        Results{c}{2,2} = cellfun(@round,Mean_Brightness_Single_RNA{c});
        end
        if exist('ci2_all')
        Results{c}{2,3} = cellfun(@round,ci2_all{c}(1,1),'UniformOutput',false);
        end
        if exist('Number_of_RNAS_per_Condensate')
        Results{c}{2,4} = Number_of_RNAS_per_Condensate{c};
        end
        if exist('Std_Inf_Number_of_RNAS_per_Condensate')
        Results{c}{2,5} = Std_Inf_Number_of_RNAS_per_Condensate{c};
        end
        if exist('Std_Sup_Number_of_RNAS_per_Condensate')
        Results{c}{2,6} = Std_Sup_Number_of_RNAS_per_Condensate{c};
        end
        if exist('Condensate_Region_Size')
        Results{c}{2,7} = Condensate_Region_Size{c};
        end
        if exist('Condensate_Int_Region_Size')
        Results{c}{2,8} = Condensate_Int_Region_Size{c};
        end
        if exist('RNAs_Density_in_Condensates')
        Results{c}{2,9} =  RNAs_Density_in_Condensates{c};
        end
        if exist('Percentage_of_RNAS_inside_Condensates_Collapse_images')
        Results{c}{2,10} = Percentage_of_RNAS_inside_Condensates_Collapse_images{c};
        end
        if exist('Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_images')
        Results{c}{2,11} = Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_images{c};
        end
        if exist('Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_images')
        Results{c}{2,12} = Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_images{c};
        end
        if exist('Localization_Normalized_Intensity_RNAs')
        Results{c}{2,13} = Localization_Normalized_Intensity_RNAs{c};
        end
        if exist('Localization_Normalized_Intensity_Condensates')
        Results{c}{2,14} = Localization_Normalized_Intensity_Condensates{c};
        end
        if exist('Number_of_Localizations')
        Results{c}{2,15} = Number_of_Localizations(c);
        end
        if exist('Number_of_Localizations_perStack')
        Results{c}{2,16} = Number_of_Localizations_perStack{c};
        end

      

        mkdir(strcat(startlocn{1},'\Results'));
        save(fullfile(strcat(startlocn{1},'\Results'),strcat('RNA_Counting_Results','.mat')),'Results');
      end

        Results_Ensemble{1,1} = 'Condition';
        Results_Ensemble{1,2} = 'Percentage of RNA Inside Condensates (0 to 1) - Ensemble (All images and Z-stacks)';
        Results_Ensemble{1,3} = 'Percentage of RNA Inside Condensates Inferior 95 Confidence Interval (0 to 1) - Ensemble (All images and Z-stacks)';
        Results_Ensemble{1,4} = 'Percentage of RNA Inside Condensates Superior 95 Confidence Interval (0 to 1) - Ensemble (All images and Z-stacks)';
        Results_Ensemble{1,5} = 'Normalized RNAs Intensities (Outside Condensates) - Ensemble (All images and Z-stacks)';     
        Results_Ensemble{1,6} = 'Normalized RNAs Intensities (Inside Condensates) (Only Positive Condensates) - Ensemble (All images and Z-stacks)'; 
        Results_Ensemble{1,7} = 'RNAs Intensities (Outside Condensates) - Ensemble (All images and Z-stacks)';     
        Results_Ensemble{1,8} = 'RNAs Intensities (Inside Condensates) (Only Positive Condensates) - Ensemble (All images and Z-stacks)'; 
        
        if exist('Condition')
        Results_Ensemble{2,1} = Condition;
        end
        if exist('Percentage_of_RNAS_inside_Condensates_Collapse_all')
        Results_Ensemble{2,2} = Percentage_of_RNAS_inside_Condensates_Collapse_all;
        end
        if exist('Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_all')
        Results_Ensemble{2,3} = Std_Inf_Percentage_of_RNAS_inside_Condensates_Collapse_all;
        end
        if exist('Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_all')
        Results_Ensemble{2,4} = Std_Sup_Percentage_of_RNAS_inside_Condensates_Collapse_all;
        end
        if exist('Localization_Normalized_Intensity_RNAs_Collapse_all')
        Results_Ensemble{2,5} = Localization_Normalized_Intensity_RNAs_Collapse_all;
        end
        if exist('Localization_Normalized_Intensity_Condensates_Collapse_all')
        Results_Ensemble{2,6} = Localization_Normalized_Intensity_Condensates_Collapse_all;
        end
        if exist('Localization_Intensity_RNAs_Collapse_all')
        Results_Ensemble{2,7} = Localization_Intensity_RNAs_Collapse_all;
        end
        if exist('Localization_Intensity_Condensates_Collapse_all')
        Results_Ensemble{2,8} = Localization_Intensity_Condensates_Collapse_all;
        end

        save(fullfile(strcat(startlocn{1},'\Results'),strcat('RNA_Counting_Results_Ensemble','.mat')),'Results_Ensemble');
      
end


if Save_Figures == 1;
mkdir(strcat(startlocn{1},'\Results_Figures_RAW'));
FolderName = strcat(startlocn{1},'\Results_Figures_RAW');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList);
  FigHandle = FigList(iFig);
  FigName   = strcat('Figure_',num2str(iFig));
  savefig(FigHandle, fullfile(FolderName, [FigName, '.fig']));
  saveas(FigHandle, fullfile(FolderName, [FigName, '.tif']));

end

end







