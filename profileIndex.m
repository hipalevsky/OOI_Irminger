function [profile_index,updown_index] = profileIndex(depth_in)
%Input depth_in is series of depths where measurements were made. A profile
%is defined as a series of depth points constantly increasing or
%decreasing. Profile breaks occur when the the depths switch between
%increasing/decreasing.

 ddiff = diff(depth_in);
 profile_index = NaN*ones(length(depth_in),1); %initialize vector
 updown_index = NaN*ones(1, length(depth_in)); %initialize vector
    profile_index(1:2) = 1; %begin with profile_index = 1;
    
 for i = 3:length(depth_in)
     %Check if ddiff(i-1), which is index for depth_in(i) and depth_in(i-1) comparison has same sign as ddiff(i-2)
     if ddiff(i-1) < 0 & ddiff(i-2) < 0 | ddiff(i-1) > 0 & ddiff(i-2) > 0
         profile_index(i) = profile_index(i-1);
     else
         profile_index(i) = profile_index(i-1) + 1;
     end
     if ddiff(i-1) < 0
         updown_index(i) = 1;
     elseif ddiff(i-1) > 0
         updown_index(i) = -1;
     end
 end
 
end