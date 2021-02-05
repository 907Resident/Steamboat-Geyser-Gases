
%% NGBN_STMB_DataAnalysis

% NGBN_STMB_DataAnalysis
% Geoscientists: Ajayi, M. & Mans, W.
% Location: Norris Geyser Basin - Steamboat (NGBN), Yellowstone National
% Park, WY
% GPS Coordinates: [44.723408, -110.703405]
% Data Type: Transect Measurements 
% Metadata File: "...\Yellowstone\Jun2019\12Jun2019\NGBN_STMB_Metadata"

%% Extra Functions Used in this Script
    % The following functions that were made by Ajayi, M.
    
        % lin_regress()
        % KeelingCurve()
        % StatsPlot_lite()
        % CH4_vs_Time()
        % CO2_vs_Time()
        % fullscreen()
        
%% Import Data Pre-Proccessed Data

% Data should be imported via "Merge_PicData_Fxn"
 directory = "E:\moyoa\Documents\OneDrive - Vanderbilt\PhD_Dissertation\Data_Analysis\Picarro\Yellowstone\Jun2019\12Jun2019\*.dat";
[TT_PicData,PD_mtrx_nm,PD_mtrx_dt, dte] =                               ...
 Merge_PicData_Fxn(directory, 'America/Denver');  

%% Gather the NGBN_STMB Data from the Eosense Matrix 'PD_mtrx'

% Calculate the relative time for each measurement 
    % Pre-allocate vector for time normalized in milliseconds.
    % Milliseconds are used because of greater precision than seconds alone
    ms  = zeros([length(TT_PicData.Time), 1]);
    % Use for-loop to calculate the difference in times between each
    % measurement
    for i = 2:length(TT_PicData.Time)
            ms(i)               = ms(i-1) +                             ...
            ( milliseconds(TT_PicData.Time(i) - TT_PicData.Time(i-1)) );
    end
    % Convert milliseconds back to seconds
    s                           = ms ./ 1E03;
    NGBN_STMB_rel_sec      = s;
    % Add the relative seconds vector to the timetable
    TT_PicData          = addvars(TT_PicData,NGBN_STMB_rel_sec,    ...
                        'Before','HP_CH4_dry');
    % Clear extraneous variables
    clearvars ms s

%% Spec Check
    % Ensure that measurements are not outside of the Picarro specs (e.g.
    % [CH4] <= 10 ppm & [CO2] <= 2500 ppm)

% Filter out Erroneous Measurements
TT_PicData.HP_CH4_dry(TT_PicData.HP_CH4_dry >= 1.00E01)         = NaN;
TT_PicData.HR_CH4_dry(TT_PicData.HR_CH4_dry >= 1.00E01)         = NaN;
TT_PicData.HP_d13_CH4(TT_PicData.HP_d13_CH4 >= 1.00E06)         = NaN;
TT_PicData.CO2_dry(TT_PicData.CO2_dry       >= 2.50E03)         = NaN;
TT_PicData.d13_CO2(TT_PicData.d13_CO2       >= 1.00E06)         = NaN;                             

%% Incorporate User Options

writeXL = input('Enter ''1'' to write main variables to an Excel file, enter zero if no Excel files are to be wrtitten: ');
% State how many times the eosAC cycled (opened and closed) during the
% long-term measurement
  % This will help to separate each closed chamber period during
  % that measurement session
nchams  = input('Enter the number of points measured at this site (reponses must be greater than zero): ');
    if nchams <= 0
       error('Please select a value greater than zero')
    elseif mod(nchams,1) ~= 0
        error('Please select an interger value')
    end
    
%%  Get Rid of the Duplicate Measurement Values
    % In some instances, users of the Picarro have observed duplicate
    % measurements. That is, measurements over 1-3 seconds that are exactly
    % the same.  Thus, it becomes necessary to remove these duplicate
    % measurements.

% The removal of the duplicate measurements can be accomplished by creating
% a vector 1's and 0's.  The 1's will be included in rows that are kept and
% the 0's will indicate rows that are to be removed.
    % For the NGBN_STMB, the data presented groups of four thus, only
    % every fourth row is retainted
RN = zeros([length(TT_PicData.HP_CH4_dry), 1]);
    % Assigns '1' to every fourth row in the vector RN
for i = 1:4:length(RN)
RN(i) = 1;
end
    % Create a new matrix that will contain the desired measurements and
    % the keep/delete indicator by horizontally concatenating RN to the
    % measurement vectors
RN_CH4 = [RN, datenum(TT_PicData.Time), NGBN_STMB_rel_sec          ...
         TT_PicData.HP_CH4_dry,         TT_PicData.HR_CH4_dry,          ...
         TT_PicData.HP_d13_CH4,         TT_PicData.HR_d13_CH4           ... 
         TT_PicData.CO2_dry,            TT_PicData.d13_CO2,             ...
         TT_PicData.Alarm_Status];
    % Retain every fourth entry and replace other entries with NaNs
for i = 1:length(RN)
    if RN_CH4(i,1) == 0
       RN_CH4(i,:) = NaN;
    end
end
    % Delete the rows that contain NaNs (as determined by the indicator
    % vector)
RN_CH4(any(isnan(RN_CH4),2),:)  = [];
    % Create a new vector that converts the edited time serial numbers to
    % datetime objects
Time                            = datestr(RN_CH4(:,2), -1);
Time                            = datetime(Time);

% Reassign the variables so that they represent the filtered information
    % Set up variables names
    VariableNames = {'Rel_Sec', 'HP_CH4_dry', 'HR_CH4_dry', 'HP_d13_CH4',  ...
                     'HR_d13_CH4', 'CO2_dry', 'd13_CO2', 'Alarm_Status'};
    % Recreate Timetable with pared down data, now called
    % *TT_NGBN_STMB*
    TT_NGBN_STMB = timetable(Time, RN_CH4(:,03),  RN_CH4(:,04),    ...
                           RN_CH4(:,05), RN_CH4(:,06), RN_CH4(:,07),    ...
                           RN_CH4(:,08), RN_CH4(:,09), RN_CH4(:,10),    ...
                          'VariableNames', VariableNames);
    % Record the summary of the data (includes descriptive statistics)
    smry_NGBN_STMB = summary(TT_NGBN_STMB);
    % Clear extraneous variables
    clearvars RN RN_CH4 RN_Time

%------Figure 01-------%
fullscreen,    
%------Figure 1A-------%
subplot(1,2,1)
    plot(TT_NGBN_STMB.Time, TT_NGBN_STMB.HP_CH4_dry, 'o')
        title('NGBN STMB 12Jun2019 [CH_4] vs. Timestamp (ALL)')
        ylabel('[CH_4] (ppm)')
        grid on
%         legend('[CH_4]', '± 1 std', 'Mean', 'Max', 'Min', ...
%            'Location', 'Northwest')

%------Figure 1B-------%
subplot(1,2,2)
    plot(TT_NGBN_STMB.Time, TT_NGBN_STMB.CO2_dry, 'mo')
        title('NGBN STMB 12Jun2019 [CO_2] vs. Timestamp (ALL)')
        ylabel('[CO_2] (ppm)')
        grid on
%     legend('[CO_2]', '± 1 std', 'Mean', 'Max', 'Min', ...
%            'Location', 'Northwest')

%------Figure 02-------%
fig = figure;
left_ax_color          = [0 0.477 0.741];
right_ax_color         = [1 0 1];

set(fig,'defaultAxesColorOrder',[left_ax_color; right_ax_color]);

yyaxis left
% --- CH4 --- %
   h1 = plot(TT_NGBN_STMB.Time,TT_NGBN_STMB.HP_CH4_dry, 'o');
                    title('NGBN STMB 12Jun2019 [CH_4] and [CO_2] vs. Time')
                    grid on
                    % Modify the line markers
                    h1.Marker               = 'o';
                    h1.LineStyle            = 'none';
                    ylabel('[CH_4] (ppm)')
yyaxis right                    
% --- CO2 --- %
   h2 = plot(TT_NGBN_STMB.Time,TT_NGBN_STMB.CO2_dry, 'o');
                    h2.Marker               = 'o';
                    h2.LineStyle            = 'none';
                    h2.Color                = [1 0 1];
                    ylabel('[CO_2] (ppm)')

                    % Add a legend
                    legend({'[CH_4]', '[CO_2]'}, 'EdgeColor', 'none',   ...
                            'Color','none', 'Location', 'Best')
                        
% Clear extraneous variables
clearvars h1 h2 left_ax_color right_ax_color
clearvars PicData PD_mtrx_dt TT_PicData 

%% Extract and Isolate Pre/ChamON/Post Data
    % Here, we separate the data collected into packages for before,
    % during, and after chamber enclosure.  This requires visual inspection
    % of the CO2 data (see first two graphs) and comparison to the field
    % notes associated with this measurement session

% Create array of time contraints 
DATE                            = dte(1:nchams);

%----PreBack----%
    % Pre-allocate the datetime array
PreBack_Times                   = NaT([nchams, 2]);

PreBack_Times(:,1)              = datetime(                             ...
['10:13:19'; '10:35:00'; '10:54:58'; '11:14:52'; '11:34:56'; '12:14:58';...
 '12:14:58'; '12:34:45'; '12:54:52'; '13:11:46'; '13:32:10'; '13:51:43';...
 '14:11:54'; '14:31:51'; '14:51:42'; '15:11:42'; '15:32:14'],           ...                                                           ...
 'InputFormat', 'HH:mm:ss');
    % Extract the time from the combined date-time value
PreBack_Times(:,1)              = DATE + timeofday(PreBack_Times(:,1));

PreBack_Times(:,2)              = datetime(                             ...
['10:18:06'; '10:39:22'; '10:59:16'; '11:18:38'; '11:38:56'; '11:58:49';...
 '12:18:56'; '12:39:20'; '12:56:09'; '13:15:47'; '13:35:57'; '13:35:57';...
 '13:55:19'; '14:15:23'; '14:35:28'; '15:15:38'; '15:35:53'],           ...
 'InputFormat', 'HH:mm:ss');
    % Extract the time from the combined date-time value
PreBack_Times(:,2)              = DATE + timeofday(PreBack_Times(:,2));

%----ChamON----%
ChamON_Times                    = NaT([nchams, 2]);

ChamON_Times(:,1)               = datetime(                             ...
['10:18:10'; '10:39:26'; '10:59:19'; '11:18:41'; '11:39:00'; '11:58:53';...
 '12:19:00'; '12:39:24'; '12:56:12'; '13:15:50'; '13:36:00'; '13:55:22';...
 '14:15:26'; '14:35:31'; '14:55:12'; '15:15:42'; '15:35:57'],           ...
 'InputFormat', 'HH:mm:ss');

ChamON_Times(:,1)               = DATE + timeofday(ChamON_Times(:,1));

% Fill in with the **"end"** time of each chamber closure period (i.e. the
% time when the chamber opens
ChamON_Times(:,2)               = datetime(                             ...
['10:34:03'; '10:54:03'; '11:14:00'; '11:34:00'; '11:53:59'; '12:13:53';...
 '12:33:54'; '12:53:57'; '13:10:47'; '13:30:56'; '13:50:47'; '14:10:51';...
 '14:30:45'; '14:50:43'; '15:10:42'; '15:30:54'; '15:50:46'],           ...
 'InputFormat', 'HH:mm:ss');

ChamON_Times(:,2)               = DATE + timeofday(ChamON_Times(:,2));

%----PostBack----%
    % No post background data needed here; the pre background is used as
    % pre and post for these continuous measurements
    
    % Clear extraneous variables
    clearvars DATE
    
%% Extract the Individual Measurement
    % Use the start and end time from the previous section, use the
    % TIMERANGE() function to isolate the data.  Then transfer the results
    % from a timetable to a multidimensional array

    NGBN_STMB_PreBack    = NaN([100, width(TT_NGBN_STMB)+1, nchams]);
    NGBN_STMB_ChamON     = NaN([420, width(TT_NGBN_STMB)+1, nchams]);
    
for i = 1:nchams
    % PreBack
    S_Pre           = timerange(PreBack_Times(i,1), PreBack_Times(i,2), ...
                               'closed');
    T               = TT_NGBN_STMB(S_Pre,:);
    A               = NaN([height(T), width(T)+1]);
    dnum            = datenum(T.Time);
    A(:,1)          = dnum;
    A(:,2:end)      = table2array(T);
   [m, n]           = size(A);
    
    NGBN_STMB_PreBack(1:m,1:n,i) = A;
    
    % For posterity clear extraneous variables here, not required but in
    % case future changes to the code are made, this line will eliminate
    % mix ups
    clearvars T A dnum m n
    
    % ChamON
    S_ChamON    = timerange(ChamON_Times(i,1), ChamON_Times(i,2),       ...
                            'closed');
    T           = TT_NGBN_STMB(S_ChamON,:);
    A           = NaN([height(T), width(T)+1]);

    dnum        = datenum(T.Time);
    A(:,1)      = dnum;
    A(:,2:end)  = table2array(T);
    % Change the seconds to relative seconds of each individual measurement
    % (0 to ~1000/1200)
        % Calculate the differences in time (sec) between each measurement
        dmy_diff= diff(A(:,2));
        % Usimg CUMSUM(), calculate the accumlated time for each
        % measurement.  In other words, this will now show the relative
        % within a measurement period (e.g. the measurement occured 45 sec
        % into the measuring period)
        A(2:end,2)  = cumsum(dmy_diff);
        % Set the first time component to zero because no time has elapsed
        % at the first measurement within the time peropd
        A(1,2)  = 0;
    % Acquire the dimnensions of the array A so that it can be properly
    % placed in the ___ data
   [m, n]       = size(A);
    
    NGBN_STMB_ChamON(1:m,1:n,i) = A;  
    
    % Clear extra variables
    clearvars T A dnum m n dmy_diff
    
end

%% Extract Datetimes

    % PreBack
        % Pre-allocate
    NGBN_STMB_PreBack_DateTimes     = NaT([100 1 nchams]);
for idx = 1:nchams
    f                               = find(isnan(NGBN_STMB_PreBack(:,1,idx)));
    NGBN_STMB_PreBack_DateTimes(1:min(f)-1,1,idx)    = datetime(        ...
         datestr(NGBN_STMB_PreBack(1:min(f)-1,1,idx)));
end

    % ChamON
        % Pre-allocate
    NGBN_STMB_ChamON_DateTimes      = NaT([350 1 nchams]);
for idx = 1:nchams
    f                               = find(isnan(NGBN_STMB_ChamON(:,1,idx)));
    NGBN_STMB_ChamON_DateTimes(1:min(f)-1,1,idx)    = datetime(         ...
         datestr(NGBN_STMB_ChamON(1:min(f)-1,1,idx)));
end

% Delete Extraneous Variables
clearvars S_Pre S_ChamON

%% Prelinary Visualizations

% Designate the number of rows and columns for the subplots 
if nchams           <= 04
   dims_row          = 1;
   dims_col          = nchmams;  
elseif nchams       >= 04 && nchams <= 08
    dims_row         = 4;
    dims_col         = ceil(nchams ./ 2);
elseif nchams       >= 09 && nchams <= 12
    dims_row         = 8;
    dims_col         = ceil(nchams ./ 2);
else
    dims_row         = 12;
    dims_col         = ceil(nchams ./ 4);
end

%-----Figure 03-----%
% Create a (2*dims_row) x (dims_col) array scatter plots of concentration
% vs time (rel sec) for both gases
    fullscreen,
for idx = 1:nchams
    
%--CH4--%
    % This if-else block helps to determine where to place each of the
    % graphs 
    if idx <= dims_col 
        CH4_pl = idx;
    else
        CH4_pl = idx+dims_col;
    end
    subplot(dims_row,dims_col,CH4_pl)
    s_CH4 = scatter(NGBN_STMB_ChamON(:,2,idx),NGBN_STMB_ChamON(:,3,idx));
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
            if idx     <= dims_col
                trans   = 3;
                pnt     = idx;
            elseif idx > dims_col 
                trans   = 2;
                pnt     = idx-dims_col;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       'NGBN %d.%d 12Jun2019[CH_4] vs. Rel Time',      ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
            % The maximum value on the x-axis will be 1200 (= 60 sec *
            % 20 min)
            xlim([0 1200])
            % Provide tick marks in 300 sec (5 min) intervals
            xticks(0:300:1200)
            % Add labels to axis
            xlabel('Relative Time (sec)', 'FontSize', 8)
            ylabel('[CH_4] (ppm)', 'FontSize', 8)
            % Apply gridding on both axis to see the data better
            grid on
            
%--CO2--%
    % This if-else block helps to determine where to place each of the
    % graphs
    if idx <= dims_col 
        CO2_pl = idx+dims_col;
    else
        CO2_pl = idx+(2*dims_col);
    end 
    subplot(dims_row,dims_col,CO2_pl)
    s_CO2 = scatter(NGBN_STMB_ChamON(:,2,idx),NGBN_STMB_ChamON(:,7,idx));
    s_CO2.MarkerEdgeColor = 'm';
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
            if idx     <= dims_col
                trans   = 3;
                pnt     = idx;
            elseif idx > dims_col
                trans   = 2;
                pnt     = idx-dims_col;
            else
                trans   = 99;
                pnt     = 99;
            end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       'NGBN %d.%d 12Jun2019[CO_2] vs. Rel Time',       ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
            % The maximum value on the x-axis will be 1200 (= 60 sec *
            % 20 min)
            xlim([0 1200])
            % Provide tick marks in 300 sec (5 min) intervals
            xticks(0:300:1200)
            % Add labels to axis
            xlabel('Relative Time (sec)', 'FontSize', 8)
            ylabel('[CO_2] (ppm)', 'FontSize', 8)
            % Apply gridding on both axis to see the data better
            grid on        
                
end

% Create a dims_row x dims_col array scatter plots of concentration vs time
% (rel sec)

% ---- Figure 04 ---- %
    % CH4
    fullscreen,
for idx = 1:nchams

    subplot(floor(dims_row / 2),dims_col,idx)
        s_CH4 = scatter(NGBN_STMB_ChamON(:,2,idx),NGBN_STMB_ChamON(:,3,idx));
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if idx     <= dims_col
                    trans   = 3;
                    pnt     = idx;
                elseif idx  > dims_col
                    trans   = 2;
                    pnt     = idx-dims_col;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str = sprintf(                                    ...
                           'NGBN %d.%d 12Jun2019[CH_4] vs. Rel Time',  ...
                            trans,pnt);
                title(title_str, 'FontSize', 8)
                % The maximum value on the x-axis will be 1200 (= 60 sec *
                % 20 min)
                xlim([0 1200])
                % Provide tick marks in 300 sec (5 min) intervals
                xticks(0:300:1200)
                % Add labels to axis
                xlabel('Relative Time (sec)', 'FontSize', 8)
                ylabel('[CH_4] (ppm)', 'FontSize', 8)
                % Apply gridding on both axis to see the data better
                grid on
end

% ---- Figure 05 ---- %
    % CO2
    fullscreen,
for idx = 1:nchams
    
    subplot(floor(dims_row / 2),dims_col,idx)
        s_CO2 = scatter(NGBN_STMB_ChamON(:,2,idx),NGBN_STMB_ChamON(:,7,idx));
        s_CO2.MarkerEdgeColor = 'm';
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if idx     <= dims_col
                    trans   = 1;
                    pnt     = idx;
                elseif idx  > dims_col
                    trans   = 2;
                    pnt     = idx-dims_col;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str = sprintf(                                    ...
                           'NGBN %d.%d 12Jun2019[CO_2] vs. Rel Time',  ...
                            trans,pnt);
                title(title_str, 'FontSize', 8)
                % The maximum value on the x-axis will be 1200 (= 60 sec *
                % 20 min)
                xlim([0 1200])
                % Provide tick marks in 300 sec (5 min) intervals
                xticks(0:300:1200)
                % Add labels to axis
                xlabel('Relative Time (sec)', 'FontSize', 8)
                ylabel('[CO_2] (ppm)', 'FontSize', 8)
                % Apply gridding on both axis to see the data better
                grid on        
end

% Clear extraneous variables
clearvars CH4_pl CO2_pl title_str pnt trans
    
%% Spec Check
    % Ensure that measurements are not outside of the Picarro specs (e.g.
    % [CH4] <= 10 ppm & [CO2] <= 2500 ppm)

% CH4 Spec Check
for sp_idx = 1:nchams
    CH4_Sp_Chk = find(NGBN_STMB_ChamON(:,3,sp_idx) > 10);
    if  any(CH4_Sp_Chk)
        NGBN_STMB_ChamON(CH4_Sp_Chk,3,sp_idx) = NaN;
        NGBN_STMB_ChamON(CH4_Sp_Chk,5,sp_idx) = NaN;
    else
        Spec_Chk_str = sprintf('There were no [CH4] measurements greater than 10 ppm for measurement number %d.', ...
                                                                   sp_idx);
        disp(Spec_Chk_str)                                                       
    end
end

% CO2 Spec Check
for sp_idx = 1:nchams
    CO2_Sp_Chk = find(NGBN_STMB_ChamON(:,7,sp_idx) > 2500);
    if  any(CO2_Sp_Chk)
        NGBN_STMB_ChamON(CO2_Sp_Chk,7:8,sp_idx) = NaN;
    else
        Spec_Chk_str = sprintf('There were no [CO2] measurements greater than 2500 ppm for measurement number %d.', ...
                                                                   sp_idx);
        disp(Spec_Chk_str)                                                       
    end
end

clearvars CH4_Sp_Chk CO2_Sp_Chk sp_idx

%% Convert gas measurements from ppm tg cgs units

% Converting constants
Vm                              = 22.71108;  % Standard Molar Volume (mol L^-1)
CH4_molec_wt                    = 16.043;    % g/mol
CO2_molec_wt                    = 44.009;    % g/mol
% Imported concentrations (resampled)
c_CH4_ppm                       = NGBN_STMB_ChamON(:,3,:);   
c_CO2_ppm                       = NGBN_STMB_ChamON(:,7,:);
% Chamber Height (m)
    % 5in = 12.70 cm = 0.127 m
H                               = 0.127; 
% Number of observations (N)
N                               = length(c_CH4_ppm);

%%% To get fluxes, think of concentration expressed in micrograms per cubic
%%% meter ug/m3
    %%% Convert from ratio of gas to air (ppm) to mass per volume (mg /
    %%% m^3)
        %%% (Molec_wt / Vm) * (1000 L / m^3) * ppm = ug / m^3
        c_CH4_cgs               = ((CH4_molec_wt / Vm) .* 1E03)         ...
                                    .* c_CH4_ppm;
        %%%% (ug / m^3) / 1000 = mg / m^3 
        c_CH4_cgs               = c_CH4_cgs ./ 1E03;
    %%% Repeat for carbon dioxide     
        %%% (Molec_wt / Vm) * (1000 L / m^3) * ppm = ug / m^3
        c_CO2_cgs               = ((CO2_molec_wt / Vm)  .* 1E03)        ...
                                    .* c_CO2_ppm;
        %%%% (ug / m^3) / 1000 =  mg / m^3 
        c_CO2_cgs               = c_CO2_cgs ./ 1E03;        

% Add the converted cgs concentrations to the main numeric array
 for idx = 1:nchams
        NGBN_STMB_ChamON(:,10,idx)    = c_CH4_cgs(:,1,idx);
        NGBN_STMB_ChamON(:,11,idx)    = c_CO2_cgs(:,1,idx);
 end

% Clear extraneous variables
clearvars N Vm CH4_molec_wt CO2_molec_wt c_CH4_cgs_UnEx c_CO2_cgs_UnEx c_CH4_ppm c_CO2_ppm 

%% Flux Calculation

%----Linear Model Analysis----% 
    % A two-minute equilibration period (approx 33 measurements) is allowed
    % before an assessment of the linear fit.  After this two-minute
    % period, the next three minutes are analyzed (approx. 52 measurements)
    % using a simple linear regression model
        % C_t = b + m*t, where the slope, m, is dC_dt
        
    % Pre-allocate vector for fluxes (per measurement cycle)
    lin_mdl     = cell ([1          nchams      2           ]); 
    lin_slope   = zeros([2          4           nchams     2]);
    lin_flux    = zeros([nchams     1           2           ]);
    
    % Measurement start time
    flx_srt_time= NaT([nchams 1]);    
    % Measurement end time
    flx_end_time= NaT([nchams 1]);
    
    % Use for-loop to run a linear regression for each of the chamber
    % measuerments
    for i = 1:nchams
       %--CH4--%
       
       % To account for the equilibration period at the beginning of the
       % closure period, we bypass the first 120 sec (2 min) of the
       % measurement
            % start
       srt_idx_CH4         = 15;
            % end
       end_idx_CH4         = find(NGBN_STMB_ChamON(:,2,i) <= 240,1, 'last');
       
       % Generate Linear Model
            % Here I chose to go from the start of the of the timeseries to
            % the last non-NaN measurement because all of these
            % measurements were so short on account of the H2O alarms and
            % when CO2 breached 2500 ppm.  No decision has been made if CH4
            % should use data that goes past the CO2 threshold yet (24 Jan
            % 2019)
       lin_mdl{1,i,1}      = fitlm(NGBN_STMB_ChamON(15:end_idx_CH4,02,i),...
                                   NGBN_STMB_ChamON(15:end_idx_CH4,10,i)); 
       % Store pertinent coeffieicents from the linear model                        
       lin_slope(:,:,i,1)  = table2array(lin_mdl{1,i,1}.Coefficients);
       % Quantify the flux for each measurement (mg m^-2 hr^-1)
            % The the "3600" represents one hour in seconds (60 sec * 60
            % min)
       lin_flux(i,1,1)     = lin_slope(2,1,i,1) .* H .* 3600;
       
       %--CO2--%
       
       % To account for the equilibration period at the beginning of the
       % closure period, we bypass the first 120 sec (2 min) of the
       % measurement  
            % start
       srt_idx_CO2         = 15;
            % end       
       end_idx_CO2         = find(NGBN_STMB_ChamON(:,2,i) <= 240,1, 'last');
       
       % Generate Linear Model
       lin_mdl{1,i,2}      = fitlm(NGBN_STMB_ChamON(15:end_idx_CO2,02,i),...
                                   NGBN_STMB_ChamON(15:end_idx_CO2,11,i)); 
       % Store pertinent coeffieicents from the linear model                        
       lin_slope(:,:,i,2)  = table2array(lin_mdl{1,i,2}.Coefficients);
       % Quantify the flux for each measurement (mg m^-2 hr^-1)
       lin_flux(i,1,2)     = lin_slope(2,1,i,2) .* H .* 3600;
       
       % Pull out the flux measurement end time
            % start
       flx_srt_time(i)   =                                              ...
              datetime(datestr(NGBN_STMB_ChamON(srt_idx_CO2, 1, i)));        
            % end
       flx_end_time(i)   =                                              ...
              datetime(datestr(NGBN_STMB_ChamON(end_idx_CO2, 1, i)));      
       
    end
    
%----Exponential Model Analysis----%
    % A two-minute equilibration period (approx 33-36 measurements) is
    % allowed before an assessment of the exponential fit.  After this
    % two-minute period, the remaining portion of the measurement period is
    % evaluted using an exponential fit:
        % C_t = psi + (C_0 - psi) * exp(-kappa*t)
    % Pre-allocate vector for fluxes (per measurement cycle)
    exp_mdl     = cell([1 nchams 2]);
    % Define the options for the exponential fit (i.e. fitting method,
    % starting point for the iteration)
    exp_ft_ops  = fitoptions('StartPoint', [0 0],                       ...
                             'Method', 'NonlinearLeastSquares');
    % Create an exponential model to evaluate the flux data using fittype()
    exp_fxn     = fittype(@(psi, kappa, c_0, t) psi + (c_0 - psi)       ...
                          * exp(-kappa*t),                              ...
                          'dependent', {'C_t'}, 'independent', {'t'},   ...
                          'coefficients', {'psi', 'kappa'}, 'problem',  ...
                          'c_0', 'options', exp_ft_ops);
        % (This is equivalent to 1.9 ppm CH4 and could be adjusted to
        % indicate the concentration at time = 0)
    exp_slope   = zeros([1  3       nchams      2]);
    exp_flux    = zeros([1  nchams              2]);
    
    x_expData   = cell([1   nchams  2            ]);
    y_expData   = cell([1   nchams  2            ]);

for i = 1:nchams
%---CH4---%
% Apply an exponential model to evaluate the flux data using fit()
    % The row indicies are hard coded to ensure that no NaN values are
    % analyzed; additionally, this segment of data accounts for the
    % majority of the data after the equilibriation period
[x_exp_prepData, y_exp_prepData]    =                                   ...
                 prepareCurveData(NGBN_STMB_ChamON(:,02,i),             ...
                                  NGBN_STMB_ChamON(:,10,i));   
x_expData{:, i, 1}                  = x_exp_prepData;
y_expData{:, i, 1}                  = y_exp_prepData;
exp_mdl  {1, i, 1}                  = fit(x_expData{i},y_expData{i},    ...
                                          exp_fxn,'problem',            ...
                                          NGBN_STMB_ChamON(01,10,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
exp_slope(1, 1, i, 1)               = exp_mdl{i}.psi;
exp_slope(1, 2, i, 1)               = exp_mdl{i}.kappa;
exp_slope(1, 3, i, 1)               = exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
exp_flux (1, i, 1)  =((exp_mdl{1, i, 1}.psi - exp_mdl{1, i, 1}.c_0) .*  ...
                       exp_mdl{1, i, 1}.kappa) .* H .* 3600;
                   
%---CO2---%
[x_exp_prepData, y_exp_prepData]    =                                   ...
                 prepareCurveData(NGBN_STMB_ChamON(:,02,i),             ...
                                  NGBN_STMB_ChamON(:,11,i));   
x_expData{:, i, 2}                  = x_exp_prepData;
y_expData{:, i, 2}                  = y_exp_prepData;
exp_mdl  {1, i, 2}                  = fit(x_expData{i},y_expData{i},    ...
                                          exp_fxn,'problem',            ...
                                          NGBN_STMB_ChamON(01,11,i));
% Gather the coeffient and the c_0 values for presentation and analysis
% later
exp_slope(1, 1, i, 2)               = exp_mdl{i}.psi;
exp_slope(1, 2, i, 2)               = exp_mdl{i}.kappa;
exp_slope(1, 3, i, 2)               = exp_mdl{i}.c_0;
% Quantify the flux for each measurement (mg m^-2 hr^-1)
exp_flux (1, i, 2)  =((exp_mdl{1, i, 2}.psi - exp_mdl{1, i, 2}.c_0) .*  ...
                       exp_mdl{1, i, 2}.kappa) .* H .* 3600;
end

% Reshape Exponential Flux to Match Shape of Linear Flux matrix [m x 1]
exp_flux = reshape(exp_flux,[length(exp_flux) 1, 2]);

% Visualize the fluxes
fullscreen,
    %  CH4
        subplot(1,2,1)
            bar(1:nchams, lin_flux(:,:,1))
            if any(abs(exp_flux(:,:,1) - lin_flux(:,:,1)) <= 1E02) == 1
                hold on
                bar(1:nchams,exp_flux(:,:,1), 0.5, 'FaceColor', 'c')
            else
                disp('The absolute difference between the exponential and linear models was greater than 100 mg m^-2 hr^-1 for all measuerments of CH4. Reexamine the models.')
            end
            grid on
            title('NGBN STMB 12Jun2019 | CH_4 Fluxes')
            xlabel('Location')
            ylabel('F_{CH_4} (mg m^{-2} hr^{-1})')
            legend({'Linear Model', 'Exponential Model'}, 'Location', 'SW')
            hold off
    % CO2
        subplot(1,2,2)
            bar(1:nchams, lin_flux(:,:,2), 'FaceColor', 'm')
            if any(abs(exp_flux(:,:,2) - lin_flux(:,:,2)) <= 1E02) == 1
                hold on
                bar(1:nchams,exp_flux(:,:,2), 0.5,                      ...
                    'FaceColor', [102/255 0/255 102/255])
            else
                disp('The absolute difference between the exponential and linear models was greater than 100 mg m^-2 hr^-1 for all measurements of CO2. Reexamine the models.')
            end
            grid on
            title('NGBN STMB 12Jun2019| CO_2 Fluxes')
            xlabel('Location')
            ylabel('F_{CO_2} (mg m^{-2} hr^{-1})')

%% Keeling Plots
    % Keeling plots (originally generated by Keeling (1961(?)) to assess
    % the source of the gas

% Pre-allocate vectors for KeelingCurve Results
KP_rslts                        = cell( [1 nchams   2        ]);
KP_coeffs                       = zeros([1 2        nchams  2]);
KP_coeffs_ci                    = zeros([2 2        nchams  2]);
KP_Fits                         = cell([nchams 4]);

for idx = 1:nchams
%--CH4--%
    % Prepare data for the regression (i.e. remove NaN values, shape the
    % vectors as vertical vectors, etc.) 
[xData, yData]  = prepareCurveData((1 ./ NGBN_STMB_ChamON(:,3,idx)),    ...
                                         NGBN_STMB_ChamON(:,5,idx));
    % Set up fittype and options.
KP_ft = fittype('poly1');

    % Fit model to data.
[KP_fitresult,KP_gof]           = fit(xData, yData, KP_ft);
 KP_Fits{idx,1}                 = KP_fitresult;
 KP_Fits{idx,2}                 = KP_gof.rsquare;

    % Gather the coefficients for the fit equation
fitvals = coeffvalues(KP_fitresult);
    % Gather the confidence intervals for c
    KP_rslts{1,idx,1}           = KP_fitresult;
    KP_coeffs(1,1,idx,1)        = fitvals(1);
    KP_coeffs(1,2,idx,1)        = fitvals(2);
    KP_coeffs_ci(:,:,idx,1)     = confint(KP_fitresult, 0.95);
    
%--CO2--%
[xData, yData]  = prepareCurveData((1 ./ NGBN_STMB_ChamON(:,7,idx)),    ...
                                         NGBN_STMB_ChamON(:,8,idx));
    % Set up fittype and options.
KP_ft                           = fittype('poly1');

    % Fit model to data.
[KP_fitresult,KP_gof]           = fit(xData, yData, KP_ft);
 KP_Fits{idx,3}                 = KP_fitresult;
 KP_Fits{idx,4}                 = KP_gof.rsquare;

    % Gather the coefficients for the fit equation
fitvals                         = coeffvalues(KP_fitresult);
    % Gather the confidence intervals for c
    KP_rslts{1,idx,1}           = KP_fitresult;
    KP_coeffs(1,1,idx,2)        = fitvals(1);
    KP_coeffs(1,2,idx,2)        = fitvals(2);
    KP_coeffs_ci(:,:,idx,2)     = confint(KP_fitresult, 0.95);
end

% Visualize the Keeling Plot Intercepts
    % Add the error bars for appropriate sites
fullscreen,
    %  CH4
        subplot(1,2,1)
            dmy_Y_CH4     = KP_coeffs(:,2,:,1);
            Y             = reshape(dmy_Y_CH4,[1 nchams]);
            neg_err_CH4   = reshape(KP_coeffs_ci(1,2,:,1),[1 nchams]);
            pos_err_CH4   = reshape(KP_coeffs_ci(2,2,:,1),[1 nchams]);
            bar(Y)
            hold on
            errorbar(1:nchams,Y,abs(Y - neg_err_CH4),                   ...
                                abs(Y - pos_err_CH4),                   ...
                    'LineStyle', 'none', 'Color', 'k')
            grid on
            title('NGBN STMB 12Jun2019 | \delta^{13}C-CH_4 Signature')
            xlabel('Location')
            ylabel('\delta^{13}C-CH_4 (‰)')
                     
    % CO2
        subplot(1,2,2)
            dmy_Y_CO2     = KP_coeffs(:,2,:,2);
            Y             = reshape(dmy_Y_CO2,[1 nchams]);
            neg_err_CO2   = reshape(KP_coeffs_ci(1,2,:,2), [1 nchams]);
            pos_err_CO2   = reshape(KP_coeffs_ci(2,2,:,2), [1 nchams]);
            bar(Y, 'FaceColor', 'm')
            hold on
            errorbar(1:nchams,Y,abs(Y - neg_err_CO2),                   ...
                                abs(Y - pos_err_CO2),                   ...
                    'LineStyle', 'none', 'Color', 'k')
            grid on
            title('NGBN STMB 12Jun2019 | \delta^{13}C-CO_2 Signatures')
            xlabel('Location')
            ylabel('\delta^{13}C-CO_2 (‰)','Interpreter', 'tex')
            
    % Clear extraneous variables
    clearvars dmy_Y_CH4 dmy_Y_CO2 Y KP_fitresul KP_gof neg_err_CO2 
    clearvars pos_err_CO2 neg_err_CH4 pos_err_CH4

% Visualize Keeling Plots
fullscreen,
for idx = 1:nchams
%--CH4--%
    % This if-else block helps to determine where to place each of the
    % graphs 
    if idx    <= dims_col 
        CH4_pl = idx;
    else
        CH4_pl = idx+dims_col;
    end
    subplot(dims_row,dims_col,CH4_pl)
        plot(KP_Fits{idx,1},                                               ...
                  1./NGBN_STMB_ChamON(:,3,idx),NGBN_STMB_ChamON(:,5,idx),  ...
                  'predfunc')
        fxn        = sprintf('?^{13}C-CH_4 = %2.2f*[CH_4]^{-1} + %2.2f',...
                              KP_coeffs(1,1,idx,1),KP_coeffs(1,2,idx,1));
        str        = sprintf('R^{2} = %1.3f', KP_Fits{idx,2});
        annotation('textbox',[0 0.8 0.5 0.2],'String',fxn,              ...
                   'FitBoxtoText','on', 'FontSize', 7,                  ...
                   'FontWeight', 'bold',                                ...
                   'BackgroundColor', [0.97 0.97 0.97], 'FaceAlpha',0.5,...
                   'EdgeColor', 'none')
        annotation('textbox',[0 0.5 0.3 0.2],'String',str,              ...
                   'FitBoxtoText','on' , 'FontSize', 7)
        xlabel('[CH_4]^{-1} (ppm^{-1})','FontSize', 7)
        ylabel('\delta^{13}C-CH_4 (‰)','Interpreter', 'tex',           ...
               'FontSize', 7)
        grid on
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
                if  idx    <= dims_col
                    trans   = 3;
                    pnt     = idx;
                elseif idx  > dims_col
                    trans   = 2;
                    pnt     = idx - dims_col;
                else
                    trans   = 99;
                    pnt     = 99;
                end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       'NGBN STMB %d.%d 12Jun2019 CH_4 Keeling Plot',   ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
%--CO2--%
    % This if-else block helps to determine where to place each of the
    % graphs
    if idx <= dims_col 
        CO2_pl = idx+dims_col;
    else
        CO2_pl = idx+(2*dims_col);
    end 
    subplot(dims_row,dims_col,CO2_pl)
        p_CO2 = plot(KP_Fits{idx,3},                                    ...
                1./NGBN_STMB_ChamON(:,7,idx),NGBN_STMB_ChamON(:,8,idx),'predfunc');
        p_CO2(1).Color = 'm';
        p_CO2(2).Color = 'k';
        p_CO2(3).Color = 'k';
        fxn            = sprintf('?^{13}C-CO_2 = %2.2f*[CO_2]^{-1} + %2.2f',...
                                 KP_coeffs(1,1,idx,2),KP_coeffs(1,2,idx,2));                    
        str            = sprintf('R^{2} = %1.3f', KP_Fits{idx,4});
        annotation('textbox',[0 0.2 0.5 0.2],'String',fxn,              ...
                   'FitBoxtoText','on', 'FontSize', 7,                  ...
                   'FontWeight', 'bold',                                ...
                   'BackgroundColor', [0.97 0.97 0.97], 'FaceAlpha',0.5,...
                   'EdgeColor', 'none')
        annotation('textbox',[0 0 0.3 0.2],'String',str,                ...
                   'FitBoxtoText', 'on', 'FontSize', 7)
        xlabel('[CO_2]^{-1} (ppm^{-1})','FontSize',7)
        ylabel('\delta^{13}C-CO_2 (‰)','FontSize',7)
        grid on
    % This if-elseif-else block is used to appropriately designate the
    % numerical portion of each title in the array of plots
                if  idx    <= dims_col
                    trans   = 1;
                    pnt     = idx;
                elseif idx  > dims_col
                    trans   = 1;
                    pnt     = idx;
                else
                    trans   = 99;
                    pnt     = 99;
                end
            % Title string (to be interated so the proper graph is
            % represented)
            title_str = sprintf(                                        ...
                       'NGBN STMB %d.%d 12Jun2019 CO_2 Keeling Plot',        ...
                        trans,pnt);
            title(title_str, 'FontSize', 8)
end

% Calculate alpha
    % a = [(?13C-CO2 - ?13C-CH4)/1000 ] + 1
    Alpha           = zeros([1,nchams]);
    thous_ln_alpha  = zeros([1,nchams]);
    for i = 1:nchams
        Alpha(i)            = ((KP_coeffs(1,2,i,2) - KP_coeffs(1,2,i,1))...
                                ./ 1E03) + 1;
        thous_ln_alpha(i)   = 1E03 .* log(Alpha(i));
    end

    % Calculate alpha using the method in Forde et al. 2019
        % a = [1E03 + ?13C-CO2] / [1E03 + ?13C-CH4]
    alpha_Forde19   = zeros([1 nchams]);
    for i = 1:nchams
        alpha_Forde19(i)    = (1E03 + KP_coeffs(1,2,i,2)) ./            ...
                              (1E03 + KP_coeffs(1,2,i,1) );
    end

%% Carbon Isotope Comparisons
    % Using a chart (originally generated by Whiticar (1999) we can further
    % elucidate the source of the gas via the carbon isotopes of the gas

    % Isolate the isotopic data
    iCH4 = NGBN_STMB_ChamON(:,5,:); iCO2 = NGBN_STMB_ChamON(:,8,:);
    % Use the isocomp_array() to plot the iCO2 vs. iCH4 data
    isocomp_array(iCH4, iCO2, ceil(dims_row/2), dims_col, nchams)    
    
%% CO2/CH4 Ratios

% Create an array of CH4/CO2 vs Time for each point measured
NGBN_STMB_ChamON_Ratio           = NGBN_STMB_ChamON(:,7,:) ./ NGBN_STMB_ChamON(:,3,:);
% Record the descriptive stats for the CH4/CO2 ratios
    % ... _stats = [mean; median; range; standard deviation];
    % *** COMPLETE the block BELOW (21 Aug 2018) ***
NGBN_STMB_ChamON_Ratio_stats     = [mean(NGBN_STMB_ChamON_Ratio);                 ...
                               median(NGBN_STMB_ChamON_Ratio);               ...
                               range(NGBN_STMB_ChamON_Ratio);                ...
                               std(NGBN_STMB_ChamON_Ratio)];

    % CO2/CH4
    fullscreen,
for idx = 1:nchams

    subplot(ceil(dims_row /2),dims_col,idx)
        s_CH4 = scatter(NGBN_STMB_ChamON(:,2,idx), NGBN_STMB_ChamON_Ratio(:,idx), ...
                        36, [127/255 0/255 255/255]);
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
            if idx     <= nchams / 2
                trans   = 3;
                pnt     = idx;
            elseif idx  > nchams / 2
                trans   = 2;
                pnt     = idx - dims_col;
            else
                trans   = 99;
                pnt     = 99;
            end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str   = sprintf(                                  ...
                             'NGBN STMB %d.%d 12Jun2019[CO_2]/[CH_4] vs. Rel Time',  ...
                              trans,pnt);
                title(title_str, 'FontSize', 8)
                % The maximum value on the x-axis will be 1200 (= 60 sec *
                % 20 min)
                xlim([0 1200])
                % Provide tick marks in 300 sec (5 min) intervals
                xticks(0:300:1200)
                % Add labels to axis
                xlabel('Relative Time (sec)', 'FontSize', 10)
                ylabel('[CO_2]/[CH_4]', 'FontSize', 10)
                % Apply gridding on both axis to see the data better
                grid on
end

%Clear extraneous variables
clearvars title_str trans pnt
    
%% 3-Dimensional Analysis

% Create a plot matrix of 3-D plots using d13C-CH4 (x-axis), d13C-CO2
% (y-axis), [CO2]-[CH4] ratio (z-axis)

    % Organize the data into the vertical vectors for the x-, y-, and
    % z-axis.  This will be a matrix with with a length of  , three columns
    % representing the metrics corresponding to the three axis, and a
    % third dimension in which each "page" represents an individual
    % chamber measurement.  The layout of the "xyz" matrix is:
    % [num of obs, x/y/z axis, cham meas (1:ncham)]
    xyz = NaN([length(NGBN_STMB_ChamON), 3, nchams]);
    for i = 1:nchams
        xyz(:,1,i) = NGBN_STMB_ChamON(:,5,i);
        xyz(:,2,i) = NGBN_STMB_ChamON(:,8,i);
        xyz(:,3,i) = NGBN_STMB_ChamON_Ratio(:,1,i);
    end

% Pre-allocate the surface fits for the 3-D plots
    sf3 = cell([nchams, 3]);
    % 3-D Plots
    fullscreen,
for idx = 1:nchams
    % Extract Data for each chamber measurement
    data_3D         = xyz(:,:,idx);
    
    % Eliminate the NaNs 
    data_3D_no_NaNs = data_3D(~any(isnan(data_3D),2),:);
    
    % Create a simple polynomial surface in the X & Y directions for the
    % plot3 charts
    [sf3{idx,1}, sf3{idx,2}, sf3{idx,3}]                                ...
                    = fit([data_3D_no_NaNs(:,1), data_3D_no_NaNs(:,2)], ...
                           data_3D_no_NaNs(:,3), 'poly22');
    % Create subplot matrix
    subplot(ceil(dims_row ./ 2),dims_col,idx)
        p3_CH4              =  plot(sf3{idx},                           ...
                              [data_3D_no_NaNs(:,1), data_3D_no_NaNs(:,2)],...
                               data_3D_no_NaNs(:,3));
       % Apply gridding on both axis to see the data better
                       grid on
                       
                       set(gcf, 'color', [0.80 0.80 0.80]);
                       colormap  cool
                       alpha= 0.65;
                       view(-62.14,32.56)
        % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if idx     <= dims_col
                    trans   = 3;
                    pnt     = idx;
                elseif idx  > dims_col
                    trans   = 2;
                    pnt     = idx - dims_col;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str   = sprintf(                                  ...
                             'NGBN STMB %d.%d 12Jun2019 Conc and Isotope Data',...
                              trans,pnt);
                title(title_str, 'FontSize', 7)
                % Add labels to axis
                xlabel('\delta^{13}C-CH_{4} (‰)', 'FontSize', 8)
                ylabel('\delta^{13}C-CO_{2} (‰)', 'FontSize', 8)
                zlabel('[CO_2]/[CH_4]', 'FontSize', 8)
                
end

% Plot residuals of the surface fit
fullscreen,
for i = 1:nchams
    
     subplot(ceil(dims_row ./ 2),dims_col,i)
        hist_resid          = histogram(sf3{i,3}.residuals,             ...
                              'NumBins', 8,                            ...
                              'normalization', 'probability');
    
    % This if-elseif-else block is used to appropriately designate the
        % numerical portion of each title in the array of plots
                if i       <= dims_col
                    trans   = 3;
                    pnt     = i;
                elseif i  > dims_col
                    trans   = 2;
                    pnt     = i - dims_col;
                else
                    trans   = 99;
                    pnt     = 99;
                end
                % Title string (to be interated so the proper graph is
                % represented)
                title_str   = sprintf(                                  ...
                             'NGBN STMB %d.%d 12Jun2019 Residual of Surface Fit',...
                              trans,pnt);
                title(title_str, 'FontSize', 7)
                % Add labels to axis
                xlabel('Residual', 'FontSize', 8)
                ylabel('Probability of Observance', 'FontSize', 7)
end

%Clear extraneous variables
clearvars title_str trans pnt data_3D data_3D_no_NaNs

%% Perimeter EGM Measurement


%% Export Data
    
% Write the pre-background, chamber, and post-background data into an Excel
% file to be used by non-MATLAB users who are a part of the project

if writeXL == 1
    for i = 1:nchams
        % Create a string to label the  tab of the excel sheet
            if i               <= nchams
                trans           = 1;
                pnt             = i;
            else
                trans           = 99;
                pnt             = 99;
            end
        % Title string (to be interated so the location is by the name
        % of the tab)
            export_tabname      = sprintf('NGBN_STMB_%d.%d',trans,pnt);
        
        % Export the chamber datetime objects (must import as strings to
        % allow xlswrite() to work)
        export_ChamON_times     = string(NGBN_STMB_ChamON_DateTimes(:,1,i));
        % Export the end of the flux measurement time
        export_flx_end_time     = flx_end_time;
        % Emport Chamber Data
        export_chamON_data      = NGBN_STMB_ChamON(:,2:end,i);
        % Export regression data
        export_lin_mdl_CH4      = lin_slope(:,:,i,1);
        export_lin_mdl_CO2      = lin_slope(:,:,i,2);
        export_lin_flux_CH4     = lin_flux(i,1,1);
        export_lin_flux_CO2     = lin_flux(i,1,2);
        export_exp_mdl          = exp_slope(:,:,i);
        export_exp_flux         = exp_flux(i);
        % Export Keeling Plot (KP) data
        export_KP_coeffs_CH4    = KP_coeffs(:,:,i,1);
        export_KP_coeffs_CO2    = KP_coeffs(:,:,i,2);
        export_KP_coeffs_ci_CH4 = KP_coeffs_ci(:,:,i,1);
        export_KP_coeffs_ci_CO2 = KP_coeffs_ci(:,:,i,2);
        % Export Pre- and Post-backgound times (import as strings)
        export_PreB_Times       = string(NGBN_STMB_PreBack_DateTimes(:,1,i));
%       export_PostB_Times      = string(NGBN_STMB_PostBack_DateTimes(:,1,i));
        % Export Pre- and Post-background gas data
        export_PreB_CH4         = NGBN_STMB_PreBack (:,3,i);
        export_PreB_CO2         = NGBN_STMB_PreBack (:,7,i);
        export_PreB_iCH4        = NGBN_STMB_PreBack (:,5,i);
        export_PreB_iCO2        = NGBN_STMB_PreBack (:,8,i);
%       export_PostB_CH4        = NGBN_STMB_PostBack(:,3,i);
%       export_PostB_CO2        = NGBN_STMB_PostBack(:,7,i);
%       export_PostB_iCH4       = NGBN_STMB_PostBack(:,5,i);
%       export_PostB_iCO2       = NGBN_STMB_PostBack(:,8,i);
        XL_filename             = "E:\moyoa\Documents\OneDrive - Vanderbilt\PhD_Dissertation\Data_Analysis\Picarro\Yellowstone\Jun2019\12Jun2019\NGBN_STMB_Measurements.xlsx";
        % Write the data into excel files for each variable
        XL_time_ChamON          = xlswrite(XL_filename, export_ChamON_times,    export_tabname, ' A3');
        XL_ChamON               = xlswrite(XL_filename, export_chamON_data,     export_tabname, ' B3');
        XL_lin_mdl_CH4          = xlswrite(XL_filename, export_lin_mdl_CH4,     export_tabname, ' M3');
        XL_lin_flux_CH4         = xlswrite(XL_filename, export_lin_flux_CH4,    export_tabname, ' R3');
        XL_exp_mdl              = xlswrite(XL_filename, export_exp_mdl,         export_tabname, ' T3');
        XL_exp_flux             = xlswrite(XL_filename, export_exp_flux,        export_tabname, ' X3');
        XL_lin_mdl_CO2          = xlswrite(XL_filename, export_lin_mdl_CO2,     export_tabname, ' Z3');
        XL_lin_flux_CO2         = xlswrite(XL_filename, export_lin_flux_CO2,    export_tabname, 'AE3');
        XL_KP_coeffs_CH4        = xlswrite(XL_filename, export_KP_coeffs_CH4,   export_tabname, 'AG3');
        XL_KP_coeffs_ci_CH4     = xlswrite(XL_filename, export_KP_coeffs_ci_CH4,export_tabname, 'AG4');
        XL_KP_coeffs_CO2        = xlswrite(XL_filename, export_KP_coeffs_CO2,   export_tabname, 'AJ3');
        XL_KP_coeffs_ci_CO2     = xlswrite(XL_filename, export_KP_coeffs_ci_CO2,export_tabname, 'AJ4');
        XL_time_PreB            = xlswrite(XL_filename, export_PreB_Times,      export_tabname, 'AM3');
        XL_PreB_CH4             = xlswrite(XL_filename, export_PreB_CH4,        export_tabname, 'AN3');
        XL_PreB_CO2             = xlswrite(XL_filename, export_PreB_CO2,        export_tabname, 'AP3');
        XL_PreB_iCH4            = xlswrite(XL_filename, export_PreB_iCH4,       export_tabname, 'AO3');
        XL_PreB_iCO2            = xlswrite(XL_filename, export_PreB_iCO2,       export_tabname, 'AQ3');
        XL_flx_end_time         = xlswrite(XL_filename, export_flx_end_time,    export_tabname, 'BH3');
%       XL_time_PostB           = xlswrite(XL_filename, export_PostB_Times,     export_tabname, 'AS3');
%       XL_PostB_CH4            = xlswrite(XL_filename, export_PostB_CH4,       export_tabname, 'AT3');
%       XL_PostB_CO2            = xlswrite(XL_filename, export_PostB_CO2,       export_tabname, 'AV3');
%       XL_PostB_iCH4           = xlswrite(XL_filename, export_PostB_iCH4,      export_tabname, 'AU3');
%       XL_PostB_iCO2           = xlswrite(XL_filename, export_PostB_iCO2,      export_tabname, 'AW3');
        % Export alpha value (alpha:d13C-CO2--d13C-CH4) for Horita
        % Geothermometer
        export_alpha            = Alpha(i);
        export_alpha_Forde19    = alpha_Forde19(i);  
        % Message to indicate one tab has been filled in with data
        msg_str = sprintf('Tab ''NGBN_STMB_%d'' has been created', i);
        disp(msg_str)
    end
        disp('All data has been exported to the worksheet')
end

clearvars export_tab_name         export_ChamON_times    export_chamON_data
clearvars export_lin_flux_CH4     export_lin_mdl_CO2     export_lin_flux_CH4
clearvars export_lin_flux_CO2     export_exp_mdl         export_exp_flux
clearvars export_KP_coeffs_CH4    export_KP_coeffs_CO2   export_KP_coeffs_ci_CH4
clearvars export_KP_coeffs_ci_CO2 export_PreB_Times      export_PostB_Times
clearvars export_PreB_CH4         export_PreB_CO2        export_PreB_iCH4
clearvars export_PreB_iCO2       %export_PostB_CH4       export_PostB_CO2
%clearvars export_PostB_iCH4      export_PostB_iCO2

clearvars XL_time_ChamON          XL_ChamON              XL_lin_mdl_CH4
clearvars XL_lin_flux_CH4         XL_lin_flux_CH4        XL_exp_mdl
clearvars XL_exp_flux             XL_lin_mdl_CO2         XL_lin_flux_CO2
clearvars XL_KP_coeffs_CH4        XL_KP_coeffs_ci_CH4    XL_KP_coeffs_CO21
clearvars XL_KP_coeffs_ci_CO2     XL_time_PreB           XL_PreB_CH4
clearvars XL_PreB_CO2             XL_PreB_iCH4           XL_PreB_iCO2
clearvars XL_time_PostB           XL_PostB_CH4           XL_PostB_CO2
%clearvars XL_PostB_iCH4          XL_PostB_iCO2          filename
clearvars XL_alpha                XL_alpha_Forde19       filename

%% Export a MATFILE
    % This will allow data to be manipulated in other scripts or elsewhere
    % in the MATLAB environment
    filename = "E:\moyoa\Documents\OneDrive - Vanderbilt\PhD_Dissertation\Data_Analysis\Picarro\Yellowstone\Jun2019\12Jun2019\NGBN_STMB_MATFILE.mat";
    save(filename,                                                         ...
                     'NGBN_STMB_PreBack_DateTimes' , 'NGBN_STMB_PreBack' , ...
                     'NGBN_STMB_ChamON_DateTimes'  , 'NGBN_STMB_ChamON'  , ...
                     'TT_NGBN_STMB'                , 'flx_end_time'      , ...
                     '-v7.3'); 

                        