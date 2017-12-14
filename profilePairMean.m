function [scivars_pair,ind_pair] = profilePairMean(gridded_in,tol);

%%% Function to take alternating up and down profiles and pair each up with
%%% a corresponding down profile, taking the mean

%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%
% gridded_in: structure of data already evenly gridded, containing:
    % scivars
    % time_start
    % updown
% tol: tolerance for time_start difference between profiles (days)

%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%
% scivars_pair: mean of paired up and down profiles
% ind_pair: profile index of the initial up index for each pair

%Use up followed by down profiles so that time between measurements at
%surface is smaller (since expect more variability at surface)
indup = find(gridded_in.updown > 0);

% [m,n,p] = size(gridded_in.scivars);
% scivars_pair = NaN*ones(m,n,length(indup));
j = 0; %initialize index grid

for i = 1:length(indup) %for each up profile
    %if subsequent profile after is down and time between profile starts is less than given tolerance
    if gridded_in.updown(indup(i) + 1) < 0
        if diff(gridded_in.time_start(indup(i): indup(i)+1)) < tol
            j = j + 1; %advance index grid
            scivars_pair(:,:,j) = mean(gridded_in.scivars(:,:,indup(i): indup(i)+1),3);
            ind_pair(j) = indup(i);
        end
    end
end

end