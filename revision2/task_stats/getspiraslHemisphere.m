function spiral = getspiraslHemisphere(T_session,spiral_session,label)
%% left stim
trialN1 = 0;
spiral_Lstim_contra = zeros(141,1);
spiral_Lstim_ipsi = zeros(141,1);
clear index spiral_cell_temp spiral_cell_temp2
index = (T_session.label == label & T_session.left_contrast > 0);
spiral_cell_temp = spiral_session(index,:);
if not(isempty(spiral_cell_temp))
    trialN1 = sum(index);
    spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
    for n = 1:141
        clear spiral_temp indx1 indx2
        spiral_temp =  spiral_cell_temp2{n};
        % for left stim, right hemi is contralateral side
        indx1 = (spiral_temp(:,1) > 570 );
        if not(isempty(indx1))
            spiral_Lstim_contra(n,1) = sum(indx1);
        else
            spiral_Lstim_contra(n,1) = 0;
        end
        % for left stim, left hemi is ipsilateral side
        indx2 = (spiral_temp(:,1) < 570);
        if not(isempty(indx2))
            spiral_Lstim_ipsi(n,1) = sum(indx2);
        else 
            spiral_Lstim_ipsi(n,1) = 0;
        end
    end
end
%% right stim
clear index spiral_cell_temp spiral_cell_temp2
trialN2 = 0;
spiral_Rstim_contra = zeros(141,1);
spiral_Rstim_ipsi = zeros(141,1);
index = (T_session.label == label & T_session.right_contrast > 0);
spiral_cell_temp = spiral_session(index ,:);
if not(isempty(spiral_cell_temp))
    trialN2 = sum(index);
    spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
    for n = 1:141
        clear spiral_temp indx1 indx2
        spiral_temp =  spiral_cell_temp2{n};
        % for right stim,left hemi is contralateral side
        indx1 = (spiral_temp(:,1) < 570);
        if not(isempty(indx1))
            spiral_Rstim_contra(n,1) = sum(indx1);
        else
            spiral_Rstim_contra(n,1) = 0 ;
        end
        % for right stim,right hemi is ipsilateral side
        indx2 = (spiral_temp(:,1) > 570);
        if not(isempty(indx2))
            spiral_Rstim_ipsi(n,1) = sum(indx2);
        else
            spiral_Rstim_ipsi(n,1) = 0;
        end
    end
end
%%
trialN = trialN1+trialN2;
if trialN > 0
    spiral_contra = (spiral_Lstim_contra+spiral_Rstim_contra)./(trialN1+trialN2)*35;
    spiral_ipsi = (spiral_Lstim_ipsi+spiral_Rstim_ipsi)./(trialN1+trialN2)*35;
    spiral = [spiral_contra,spiral_ipsi];
else
    spiral = zeros(141,2);
end
