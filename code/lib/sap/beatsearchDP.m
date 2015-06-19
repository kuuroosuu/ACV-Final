function Backtrack_beats = beatsearchDP(S_onset_str,indx,Backtrack_Tempo,sgsrate,tt,cand_num,search_boundary,alpha,mode)

% ===================================== Usage =====================================
% Backtrack_beats = beatsearchDP(S_onset_str,indx,Backtrack_Tempo,sgsrate,tt,cand_num,search_boundary,alpha,mode)
%  S_onset_str : The novelty curve / onset strength waveform.
%  indx : The starting index in novelty curve for searching beats.
%  Backtrack_Tempo : The tempo curve computed from previous DP.
%  sgsrate : The sampling rate of novelty curve.
%  tt : The time index of tempo curve.
%  cand_num : It indicates how many candidates will be pick from each search 
%             region.
%  search_boundary : It defines the boundaries of search region. The format 
%                    is [x y], which means that the search region is from 
%                    (x * beat period) to (y * beat period) into search direction.
%  alpha : Which controls how tightly the estimated tempo is enforced within 
%          the search range.
%  mode : Which controls the direction of beat tracking from indx.
%         (1 : Back tracking / 2 : Forward tracking)
% =================================================================================

Stage = 2;
[sr_indx,period,search_range,tracking_end_flag] = get_SearchRange(S_onset_str,indx,sgsrate,tt,Backtrack_Tempo,search_boundary,mode);
minimum_period = (60/max(Backtrack_Tempo))*sgsrate;
% Initialization
switch mode
    case 1
        Backtrack_table = zeros(cand_num,round(indx/minimum_period));
        Backtrack_table(:,1) = [1:cand_num]';
        Backtrack_ref_indx = zeros(cand_num,round(indx/minimum_period));
        Backtrack_ref_indx(:,1) = indx;
    case 2
        Remain = length(S_onset_str) - indx + 1;
        Backtrack_table = zeros(cand_num,round(Remain/minimum_period));
        Backtrack_table(:,1) = [1:cand_num]';
        Backtrack_ref_indx = zeros(cand_num,round(Remain/minimum_period));
        Backtrack_ref_indx(:,1) = indx;
end

[state1_cum_prob,state1_indx,skip] = cand_SearchRange(S_onset_str,sr_indx,period,search_range,cand_num,alpha,tracking_end_flag,mode);
if skip == true
    Backtrack_beats = [];
    return;
end
track_next_state_holder = state1_indx(find(state1_cum_prob == max(state1_cum_prob)));

while true
    [sr_indx,period,search_range,tracking_end_flag] = get_SearchRange(S_onset_str,track_next_state_holder,sgsrate,tt,Backtrack_Tempo,search_boundary,mode);
    [state2_cum_prob,state2_indx,skip] = cand_SearchRange(S_onset_str,sr_indx,period,search_range,cand_num,alpha,tracking_end_flag,mode);
%     if skip == true
%         break;
%     end
%   moved below by Frank for correctcing end point
    track_next_state_holder = state2_indx(find(state2_cum_prob == max(state2_cum_prob)));
    state2_cum_prob = repmat(state2_cum_prob,1,cand_num);
    state2_cum_prob = state2_cum_prob';
    En = repmat(state1_cum_prob,1,cand_num) + state2_cum_prob;
    [R Indx] = max(En);
    state1_cum_prob = R';
    Backtrack_table(:,Stage) = Indx';
    Backtrack_ref_indx(:,Stage) = state1_indx(Indx');
    state1_indx = state2_indx;     
    Stage = Stage + 1;    
    if skip == true
        break;
    end   
end

% Back tracking
[Final_Cand,num_cand] = Final_state_check(state1_indx,Backtrack_ref_indx,Stage);
Backtrack_lossclip =  Backtrack_PreProc(Final_Cand,num_cand,mode);
Backtrack_beats = zeros(1,Stage-2);
[Start,Indicater] = max(state1_cum_prob);
P_Indicater = Backtrack_table(Indicater,Stage-1);
Backtrack_beats(Stage-1) = Backtrack_ref_indx(P_Indicater,Stage-1);
for n = Stage-2:-1:1
    P_Indicater = Backtrack_table(P_Indicater,n);
    Backtrack_beats(n) = Backtrack_ref_indx(P_Indicater,n);
end
Backtrack_beats = [Backtrack_beats Backtrack_lossclip'];

function [search_region,period,search_range,tracking_end_flag] = get_SearchRange(S_onset_str,indx_now,sgsrate,time_tag,tempo_vector,search_boundary,mode)

tempo = sample2tempo(indx_now,sgsrate,time_tag,tempo_vector);

period = (60/tempo)*sgsrate;
switch mode
    case 1
        next_state_central_tempo = sample2tempo(indx_now - round(period),sgsrate,time_tag,tempo_vector);
    case 2
        next_state_central_tempo = sample2tempo(indx_now + round(period),sgsrate,time_tag,tempo_vector);
    otherwise
        error('Search method seeting error (1:Backward 2:Forward)');
end

dr = abs(next_state_central_tempo - tempo);

if dr <= 10
    search_range = round(search_boundary*period);  
else
    period  = (60/max([next_state_central_tempo tempo]))*sgsrate;
    search_range = round(search_boundary*period);
end

switch mode
    case 1
        search_region = [indx_now-search_range(2),indx_now-search_range(1)];
        check_end_flag = find(search_region < 0);
    case 2
        search_region = [indx_now+search_range(1),indx_now+search_range(2)];
        check_end_flag = find(search_region > length(S_onset_str));
end


if ~isempty(check_end_flag)
    tracking_end_flag = true;
else
    tracking_end_flag = false;
    Check = localmax(S_onset_str(search_region(1):search_region(2)),0);
    if isempty(find(Check,1))
        extended_search_boundary = [0.5 4];
        search_range = round(extended_search_boundary*period);
        switch mode
            case 1
                search_region = [indx_now-search_range(2),indx_now-search_range(1)];
                check_end_flag = find(search_region <= 0);
            case 2
                search_region = [indx_now+search_range(1),indx_now+search_range(2)];
                check_end_flag = find(search_region > length(S_onset_str));
        end
        if ~isempty(check_end_flag)
            tracking_end_flag = true;
        end
    end
end

function [state_prob,state_indx,skip] = cand_SearchRange(S_onset_str,sr_indx,period,search_range,cand_num,alpha,tracking_end_flag,mode)

skip = false;
range = -search_range(2):-search_range(1);
srww = exp(-0.5*((alpha*log(range/-period)).^2)); % Search Region Weighting Window
if mode == 2
    srww = fliplr(srww);
end

if tracking_end_flag == true
    switch mode
        case 1
            zpad = min(1-sr_indx(1),(search_range(2)-search_range(1)+1));
            %search_region = abs(srww) .* [zeros(1,zpad) ,S_onset_str(sr_indx(1)+zpad:sr_indx(2))];
            search_region = abs(srww) .* [zeros(1,zpad) ,S_onset_str(1:sr_indx(2))];
        case 2
            zpad = min(sr_indx(2)-length(S_onset_str),(search_range(2)-search_range(1)+1));
            search_region = abs(srww) .* [S_onset_str(sr_indx(1):end) , zeros(1,zpad)];
    end
    %added by frank for forward ending
    if zpad >  ((search_range(2)-search_range(1)+1)*3/5)
       peak_indx=[];
    else    
      peak = localmax(search_region,0);
      peak_indx = find(peak == 1);
    end 
    if isempty(peak_indx)
        state_prob = zeros(2,1);
        state_indx = zeros(2,1);
        skip = true;
        return;
    end
else
    search_region = abs(srww) .* S_onset_str(sr_indx(1):sr_indx(2));
    peak = localmax(search_region,0);
    peak_indx = find(peak == 1);
end

offset = sr_indx(1)-1;
%added by frank for forward ending
%if ((zpad >  (search_range(2)-search_range(1)+1)/2) && (max(search_region(peak_indx))< 0.7))
%if (zpad >  (search_range(2)-search_range(1)+1)/2)
%   peak_val=zeros(size(peak_indx));
%else
  peak_val = search_region(peak_indx);
%end
peak_indx = offset + peak_indx;
peak_indx = peak_indx';
peak_prob = val_normalize(peak_val,max(search_region),min(search_region));
[state_prob,state_indx] = state_gen(peak_prob,peak_indx,cand_num);

function tempo = sample2tempo(Onset_str_indx,sgsrate,time_tag,tempo_vector)

indx_time = Onset_str_indx/sgsrate;
if indx_time < min(time_tag)
    tempo = tempo_vector(1);
else if indx_time > max(time_tag)
        tempo = tempo_vector(end);
    else
        d = abs(time_tag-indx_time);
        minpos = find(d == min(d));
        clear d;
        if numel(minpos) > 1
            %changed by frank for using highest tempo in search region
            %tempo = mean(tempo_vector(minpos));
            tempo = max(tempo_vector(minpos));
        else
            tempo = tempo_vector(minpos);
        end
    end
end

function prob = val_normalize(val,hmax,hmin)

ns1 = zeros(numel(val),1);
range = hmax-hmin;

for n = 1:numel(val)
    ns1(n) = 100*(val(n)-hmin)/range;
end

prob = ns1./sum(ns1);

function [state_prob,state_indx] = state_gen(peak_prob,peak_indx,cand_num)

if numel(peak_prob) > cand_num
    [prob_sorted,rearrange_indx] = sort(peak_prob,1,'descend');
    indx_sorted = peak_indx(rearrange_indx);
    state_prob = prob_sorted(1:cand_num);
    state_indx = indx_sorted(1:cand_num);
else if numel(peak_prob) < cand_num
        state_prob = zeros(cand_num,1);
        state_indx = zeros(cand_num,1);
        state_prob(1:numel(peak_prob)) = peak_prob;
        state_indx(1:numel(peak_prob)) = peak_indx;
    else
        state_prob = peak_prob;
        state_indx = peak_indx;
    end
end

function [result,num] = Final_state_check(state1_indx,Backtrack_ref_indx,Stage)

vaild = state1_indx(find(state1_indx > 0));
check = ismember(vaild,Backtrack_ref_indx(:,Stage-1));
non_equal_indx = find(check == false);
if isempty(non_equal_indx)
    result = [];
    num = 0;
else
    non_equal_val = state1_indx(non_equal_indx);
    result = sort(non_equal_val);
    num = numel(result);
end

function Backtrack_lossclip = Backtrack_PreProc(Final_Cand,num_cand,mode)

Backtrack_lossclip = [];

switch mode
    case 1
        if num_cand > 1
            Backtrack_lossclip = Final_Cand(1:end);  
        end
    case 2
        Backtrack_lossclip = Final_Cand;      
end

