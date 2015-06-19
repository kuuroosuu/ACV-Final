function Backtrack = DP(source,penalty,mode)

% ===================================== Usage =====================================
% Backtrack = DP(source,penalty,mode)
%  source : The input data for dynamic programming, supporting matrix and 
%           sparse matrix representation in MATLAB.
%  penalty : Which is the transition penalty for finding the tempo curve.
%  mode : Which corresponds to the data type of source.
%         mode = 1 : full matrix dynamic programming
%         mode = 2 : sparse matrix dynamic programming
%  Backtrack : Which returns the most likely curve over all columns of
%              source.
% =================================================================================

% Initialization
if mode == 1
    Normalized_S = source.Cand_value/max(max(source.Cand_value));
    m = size(source.Cand_value,1);
    n = size(source.Cand_value,2);
else
    Normalized_S = source/max(max(source));
    m = size(source,1);
    n = size(source,2);
end
Start_frame = 1;
State1 = Normalized_S(:,Start_frame);

if mode == 1
    Backtrack_table = zeros(m,n-1);
    Backtrack_ref_indx = zeros(m,n-1);
    Backtrack = zeros(1,n);
    for ts = 2:n
        Peak_location_prev = source.Cand_indx(:,ts-1);
        Peak_location_current = source.Cand_indx(:,ts);
        Peak_location_prev = repmat(Peak_location_prev,1,m);
        Peak_location_current = repmat(Peak_location_current,1,m);
        Distance_Matrix = penalty .* abs(Peak_location_current - Peak_location_prev');
        Str = Normalized_S(:,ts);
        Str = repmat(Str,1,m);
        Str = Str';
        En = repmat(State1,1,m) + Str;
        En = En - Distance_Matrix;
        [R Indx] = max(En);
        State1 = R';
        Backtrack_table(:,ts-1) = Indx';
        Backtrack_ref_indx(:,ts-1) = Peak_location_prev(Indx);
    end
    
    % Back tracking (mode 1)
    [Start,Indicater] = max(State1);
    Backtrack(n) = Peak_location_current(Indicater);
    Indicater = Backtrack_table(Indicater,end);
    for ts = n-2:-1:1
        Backtrack(ts+1) = Backtrack_ref_indx(Indicater,ts+1);
        Indicater = Backtrack_table(Indicater,ts);
    end
    Backtrack(1) = Backtrack_ref_indx(Indicater,1);
    
else if mode == 2
        Backtrack_table = cell(1,n-1);
        Backtrack = zeros(1,n);
        [Indx1,dummy,Str1] = find(State1);
        while (isempty(Indx1) && isempty(Str1))
            Start_frame = Start_frame + 1;
            State1 = Normalized_S(:,Start_frame);            
            Backtrack_table{Start_frame-1} = [0,0];
            [Indx1,dummy,Str1] = find(State1);
        end
        
        for ts = Start_frame+1:n
             [Indx2,dummy,Str2] = find(Normalized_S(:,ts));
             if (isempty(Indx2) && isempty(Str2))
                 Backtrack_table{ts-1} = [Backtrack_table{ts-2}(:,2),Backtrack_table{ts-2}(:,2)];
             else
                 State2 = zeros(length(Str2),1);
                 BTt = zeros(length(Str2),1);
                 for step = 1:length(Str2)
                     D = penalty * abs(Indx1 - Indx2(step));
                     En = Str1 + Str2(step) - D;
                     [State2(step),indx] = max(En); 
                     BTt(step) = Indx1(indx);
                 end
                 Backtrack_table{ts-1} = [BTt,Indx2];
                 clear Str1;
                 clear Indx1; 
                 Indx1 = Indx2;
                 Str1 = State2;
             end
        end
        
        % Back tracking (mode 2)
        [Start,Indicater] = max(Str1);
        Link = Backtrack_table{n-1};
        Backtrack(n) = Link(Indicater,2);
        P_Indicater = Link(Indicater,1);
        for ts = n-2:-1:1
            Link = Backtrack_table{ts};
            if isequal(Link,[0,0])
                Backtrack(ts+1) = P_Indicater;
            else
                Indicater = find(P_Indicater == Link(:,2)); 
                Backtrack(ts+1) = Link(Indicater,2);
                P_Indicater = Link(Indicater,1);
            end
        end
        Backtrack(1) = P_Indicater;
    end
end








