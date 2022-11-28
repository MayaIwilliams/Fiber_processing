% Fiber photometry data processing
%
% If you want to use this code, please cite our Jove paper:
%   Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry to
%   Record Neural Activity in Freely Moving Animal. J. Vis. Exp.
%   (152), e60278, doi:10.3791/60278 (2019).

% Run section by section (Ctrl+Enter)
clear all
clc
close all



DATADIR = 'C:\Users\degoulet.m\Documents\MATLAB\Photometry_data\Cohort2' ;
ID_S = 'R59' ; % 'rxx'  
ID_F = 'FPdataR59_Left2022-09-16T15_39_08.csv' ;
events = [] ; % put your events inside this variable (workspace > double click > copy paste)
R0 = find(events(:,1) == 0) ;
events(R0, :) = [] ;
events(:, 1) = events(:,1) / 100 ;
% for iF = 1 : length(ID)
    
    %% Your data
    df = readtable([DATADIR, filesep, ID_S, filesep, ID_F]); % Change to your file
    head(df,5)

    
    % Change the next two lines depending on your data frame
    if df.LedState(1) == 0
        
        df(1, :) = [] ;

    end
    
    if df.LedState(1) == df.LedState(end)
        
        df(end, :) = [] ;   
        
    end
    
    %% Reference
    R = find(df.LedState == 1);
    raw_reference = df.Region0G(R)';
    
    %% Time 
   zdf.Timestamp = df.Timestamp - df.Timestamp(1) ;
   Ref_time = zdf.Timestamp(R)';

    %% Signal
    P =  find(df.LedState == 2);
    raw_signal = df.Region0G(P)';
    
    %% Event
    R  = find(df.Input0 == 1) ;
    event_time = zdf.Timestamp(R(1)) ;

    %% Filtering
    R1 = find(isnan(raw_reference)) ; % Detecting NaN values
    R2 = find(isnan(raw_signal)) ; % Detecting NaN values
    raw_reference([R1 R2]) = [] ; % Removing NaN values
    raw_signal([R1 R2]) = [] ; % Removing NaN values
   Ref_time([R1 R2]) = [] ; % Removing NaN values

    %% Analysis step by step
    % Smooth data
    smooth_win = 250;
    smooth_reference = movmean(raw_reference,smooth_win);
    smooth_signal = movmean(raw_signal,smooth_win);
    
    
    %% Remove slope using airPLS algorithm (airPLS.m)
    lambda = 5e9;
    order = 2;
    wep = 0.1;
    p = 0.5;
    itermax = 50;
    [reference,base_r]= airPLS(smooth_reference,lambda,order,wep,p,itermax);
    [signal,base_s]= airPLS(smooth_signal,lambda,order,wep,p,itermax);
    
    %% Remove the begining of recordings
    
    remove = 1200;
    reference = reference(remove:end);
    signal = signal(remove:end);
    Ref_time = Ref_time(remove:end);
    
    %% Standardize signals
    z_reference = (reference - median(reference)) / std(reference);
    z_signal = (signal - median(signal)) / std(signal);
    
    %% Fit reference signal to calcium signal
    % using non negative robust linear regression
    fitdata = fit(z_reference',z_signal',fittype('poly1'),'Robust','on');
    
    %% Align reference to signal
    z_reference = fitdata(z_reference)';

    %% Calculate z-score dF/F
    zdFF = z_signal - z_reference;
    
    % Plot z-score dF/F
    figure
    plot(Ref_time, zdFF,'k')
    
    color = {'c' 'b' 'k' 'g' 'r'} ;
    
    for iE = 1 : 5
    
        R = find(events(:, 2) == iE) ;
        hold on
        plot(event_time + events(R, 1), repmat(iE, [length(event_time + events(R, 1)/60) 1]), ['.', color{iE}])
        eS(1, iE) = mean(zdFF(R)) ;
        eS(2, iE) = std(zdFF(R)) ;
    end

    %save([DATADIR, filesep, ID_S, filesep,'stats.mat'], 'eS')
    %plot TTLs

%     hold on
%     RR = find(df.Input0(R(200 : end)) == 0) ;
%     TT = zdf.Timestamp(R(200 : end)) ;
%     TT(RR) = NaN ;
%     plot(TT,repmat(2, [1,length(TT)]), '.r')
 
    end
%     snr = mean(zdFF)/std(zdFF)
%     
%     
%     %saveas(gcf,['zdFF_',ID{iF},'.png'])
%     savefig(gcf, ['STN_zdFF_',ID{iF},'.fig'])
%     
%     %% Contact us
%     % If you have any questions please contact us: ekaterina.martianova.1@ulaval.ca
%     
% end

% R  = find(df.Input0 == 1) ;
% event_time = df.Timestamp(R(1)) ;
% 
% for iE = 1 : 5
%     
%     R = find(event(:, 2) == iE) ;
%     hold on
%     plot((event_time + event(R, 1))/1000, -iE - 1 : 2/length(R) : -iE + 1 - 2/length(R) , '.r')
%     
% end
