function zIdPhi = zIdPhi(rat)
%{
Implementation of zIdPhi (z-scored integrated absolute change 
in angular velocity of the head) algorithm for Vicarious Trial 
and Error (VTE) quantification proposed by Andrew E. Papale & 
Jeffrey J. Stott & Nathaniel J. Powell & Paul S. Regier & A. 
David Redish and described in "Interactions between deliberation 
and delay-discounting in rats".

This wrapper function utilizes the foaw_diff.m function adapted 
from Dr Papale Github. Original source code can be found at
https://github.com/andrewpapale/odorTrails

This wrapper function can be modified to fit your data structure 
format.
%}

% metadata structure contains subject metadata
% i.e. exper. events, maze boundaries, xy coord 
metadata = rat.metadata;

N = length(rat.trials);
vals =[N]; %

% iterate through all trials
% start of the algorithm
for index=1:N
    
    %STEP 1
    trial = rat.trials(index); % obtain trial with index
    
    % MISC: adjust xy maze coords to standard scale.
    % This was needed to our specific dataset, as each
    % subject had a different scale of xy coordinates
    adjust = trial.get('xy_adjust'); 
    
    % crop trajectory at decision point of maze
    % MISC: in our internal data structure we had
    % find_event function, that gets the time points
    % of specific events.
    entr_decision_zone = find_event(2, trial,metadata,0);
    arm_entry = find_event(3, trial, metadata,0);
    
    % check if the slice contains the data.
    % in the case of 'omission' or 'return'
    % the data slice was empty for some trials
    if(isnan(arm_entry) || isnan(entr_decision_zone))
        continue;
    end
    
    % adjust trajectory to bring data to a standard scale
    % see above for more details on this
    [x,y,~] = get_trajectory(trial,[entr_decision_zone arm_entry],adjust);
    
    
    %STEP 2
    %{
    The values of m and d are suggested by the authors
    of the algotithm
    %}
    sr = 0.02; % sampling rate
    d = 0.05; % position noise boundary
    m = 20; 
    
    dx = foaw_diff(x, sr, m, d);
    dy = foaw_diff(y, sr, m, d);
    
    %STEP 3
    % see the original paper for more details
    Phi = atan2(dx,dy);
    
    %STEP 4
    % see the original paper for more details
    Phi = unwrap(Phi);
    dPhi = foaw_diff(Phi,0.02, m, d);
    
    
    %STEP 5
    % see the original paper for more details
    IdPhi = trapz(abs(dPhi));
    
    %STEP 6 START
    % see the original paper for more details
    vals(index) = IdPhi;
    
    val = IdPhi;
    fprintf('%d: %.7f\n',index, val);
    figure
    plot(x,y);
    
%hold on;
end
% obtain z on IdPhi by z-normalizing data across all trials.
% vals array contains a single zIdPhi value per trial. Trials
% with std > 0.5 are considered VTE. Refer to the original 
% paper by Papale et.al. (2012) for more details.  
zIdPhi = zscore(vals)

% END STEP 6
% end of the algorithm
end