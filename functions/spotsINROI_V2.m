%% Extract only the tracks in the subROIs


%Locs inside are the localizations inside the subROIs. sroi_idx
%Locs outside are the localizations outside the subROIs. Outside_idx

function [sroi_idx, Outside_idx] = spotsINROI_V2(localizations,subROIs);

        Outside_idx = 1:1:size(localizations,1);
     

   if iscell(subROIs);     
        for ii=1:size(subROIs,1); %Iterate on each subROI
            
        in{ii} = inpoly(localizations,subROIs{ii});
        
        sroi_idx{ii}   = find(in{ii} == 1);
 
       
        end
        
        all_sroi_idx = vertcat(sroi_idx{:});
        Outside_idx(all_sroi_idx) = [];
     
   else
       
            
        in = inpoly(localizations,subROIs);
        
        sroi_idx{1}   = find(in == 1);
        
        Outside_idx(sroi_idx{1}) = [];  
        
   end
   
   
   
   
end  
                    