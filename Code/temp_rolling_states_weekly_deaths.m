function temp_rolling_states_weekly(run_id)

caddate1=[2020 04 20];

datenum1=datenum(caddate1);


state_id = 52; % USA
%state_id = 33; % 33- NY
%state_id = 10; % 10 - Florida

date_id = datetime(caddate1) + run_id*7;

date=datestr(date_id,'mm-dd-yy')

rolling_window_subepidemicFrameworkCoronavirusStatesUS_Deaths(state_id,date);

end

