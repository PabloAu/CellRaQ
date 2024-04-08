%% Extract only the tracks in the subROIs


%Locs inside are the localizations inside the subROIs. sroi_idx
%Locs outside are the localizations outside the subROIs. Outside_idx

function [Inside_idx, Outside_idx, subROI_idx] = spotsINROI_V3(localizations,subROIs);

        Outside_idx = 1:1:size(localizations,1);
     

   if iscell(subROIs);     
        for ii=1:size(subROIs,1); %Iterate on each subROI
            
        in{ii} = inpoly(localizations,subROIs{ii});
        
        Inside_idx{ii}   = find(in{ii} == 1);

        subROI_idx(Inside_idx{ii}) = ii;

        
        end
        
        All_Inside_idx = vertcat(Inside_idx{:});
        Outside_idx(All_Inside_idx) = [];
     
   else
       
            
        in = inpoly(localizations,subROIs);
        
        Inside_idx{1}   = find(in == 1);
        
        Outside_idx(Inside_idx{1}) = [];  
        
   end
   
   
   
   
end  
         