% calc_flgt_envelope.m
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ n_allow ] = calc_flgt_envel(naca2415,rho,text)
    load_aircraft_parameters;
    load_conversions;
    
    v_dive   = 1.5*v_cruise;                                % Dive Speed            m/s
    v_spos   = sqrt((2*wgt_max)/(naca2415.CLmax*rho*S));    % Positive Stall Speed  m/s
    v_sneg   = sqrt((2*wgt_max)/(naca2415.CLmin*rho*S));    % Negative Stall Speed  m/s
    
    % Determine velocity profile
    v_max = max([v_spos v_sneg v_dive v_cruise v_maneuver]);
    v_plot = 0:0.25:v_max;
    
    n_pos = (0.5*rho*naca2415.CLmax*v_plot.^2*S)/(wgt_max);   % Positive Load Factor Curve
    n_neg = (0.5*rho*naca2415.CLmin*v_plot.^2*S)/(wgt_max);   % Negative Load Factor Curve
    
    tmp = abs(n_pos-n1);
    [PHAA_val PHAA_ind] = min(tmp);
    tmp = abs(v_plot-v_dive);
    [PLAA_val PLAA_ind] = min(tmp);
    tmp = abs(n_neg-n2);
    [NHAA_val NHAA_ind] = min(tmp);
    tmp = abs(v_plot-v_cruise);
    [NLAA_val NLAA_ind] = min(tmp);
    
    npos_full = [n_pos(1:PHAA_ind),n1*ones(1,PLAA_ind-PHAA_ind)];
    nneg_full = [n_neg(1:NHAA_ind),n2*ones(1,NLAA_ind-NHAA_ind),linspace(n2,n3,length(v_plot)-NLAA_ind)];
    
    % Gust Load Factor
    wgt_max_eng = mass_max*kg2lbf;
    c_eng       = c*m2ft;
    S_eng       = S*(m2ft^2);
    rho_eng     = rho*(kgm3_2_slugft3);
    wing_load   = wgt_max_eng/S_eng;
    V           = [v_cruise v_dive]*ms2knots;
    Ue          = [v_gust_cruise v_gust_dive];
    
    % n Gust at Cruise/Dive Velocity
    CL_cruise = (2*wgt_max)/(rho*S*v_cruise^2);
    CL_dive   = (2*wgt_max)/(rho*S*v_dive^2);
    mu = (2*wing_load)/(rho_eng*c_eng*g_eng*naca2415.CL_alpha);
    Kg = (0.88*mu)/(5.3+mu);
    n_gust = 1 + ((Kg*V.*Ue*naca2415.CL_alpha)/(498*wing_load));
    
    ngustpos_full = [linspace(1,n_gust(1),NLAA_ind),...
                    linspace(n_gust(1),n_gust(2),length(v_plot)-NLAA_ind)];
    ngustpos_dive_full = linspace(1,n_gust(2),length(v_plot));
    ngustneg_full = [linspace(1,-n_gust(1)+2,NLAA_ind),...
                    linspace(-n_gust(1)+2,-n_gust(2)+2,length(v_plot)-NLAA_ind)];
    ngustneg_dive_full = linspace(1,-n_gust(2)+2,length(v_plot));
    
    % Calculate Allowable Envelope
    tmp = abs(n_pos - 1);   % Load Factor = 1
    [val_1 ind_1] = min(tmp);
    tmp = abs(n_neg - (-1)); % Load Factor = -1
    [val_neg1 ind_neg1] = min(tmp);
    for i = PHAA_ind:length(v_plot)
       if(ngustpos_full(i) >= npos_full(i))
           n_temp(i-PHAA_ind+1) = ngustpos_full(i);
       else
           n_temp(i-PHAA_ind+1) = npos_full(i);
       end
    end
    npos_allow = [npos_full(ind_1:PHAA_ind-1),n_temp];
    for i = NHAA_ind:length(v_plot)
       if(ngustneg_full(i) <= nneg_full(i))
           n_temp(i-NHAA_ind+1) = ngustneg_full(i);
       else
           n_temp(i-NHAA_ind+1) = nneg_full(i);
       end
    end
    nneg_allow = [nneg_full(ind_neg1:NHAA_ind-1),n_temp];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     PLOT FLIGHT ENVELOPE                        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure();
    grid on; hold on;
    % FLIGHT ENVELOPE LINES
    p1=plot(v_plot,npos_full,'b');
    plot(v_plot,nneg_full,'b');
    line([v_dive,v_dive],[n3,n1],'Color','b');
    % GUST LOAD LINES
    p2=plot(v_plot,ngustpos_full,'r');
    plot(v_plot,ngustneg_full,'r',...
         v_plot,ngustpos_dive_full,'r',v_plot,ngustneg_dive_full,'r');
    % ALLOWABLE ENVELOPE LINES
    p3=plot(v_plot(ind_1:end),npos_allow,'LineWidth',4,'Color','g');
    line([v_plot(ind_1),v_plot(ind_1)],[0,n_pos(ind_1)],'LineWidth',4,'Color','g');
    plot(v_plot(ind_neg1:end),nneg_allow,'LineWidth',4,'Color','g');
    line([v_plot(ind_neg1),v_plot(ind_neg1)],[0,n_neg(ind_neg1)],'LineWidth',4,'Color','g');
    line([v_plot(ind_1),v_plot(ind_neg1)],[0,0],'LineWidth',4,'Color','g');
    line([v_dive,v_dive],[nneg_allow(end),n1],'LineWidth',4,'Color','g');
    
    title(['Flight Envelope: ' text],'FontSize',14);
    xlabel('Velocity (m/s)','FontSize',12);   ylabel('Load Factor, n','FontSize',12);
    
    legend([p1,p2,p3],{'Maneuvering Loads','Gust Loads','Allowable Envelope'},'FontSize',12,'Location','northwest');
    print(fig,[pwd '/Flight_Envelope_Figure/' text '_Flight_Envelope'],'-djpeg');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                   EXTRACT CRITICAL LOAD FACTORS                 %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PHAA LOAD FACTOR
    if(ngustpos_full(PHAA_ind) < npos_full(PHAA_ind))
        n_allow.name(1) = cellstr('PHAA (Maneuver)');
        n_allow.n(1)    = npos_full(PHAA_ind);
        [n_allow.AoA(1) n_allow.Cd(1)] = find_n_conditions(npos_full(PHAA_ind),...
                            v_plot(PHAA_ind),rho,S,naca2415,wgt_max);
    else
        n_allow.name(1) = cellstr('PHAA (Gust)');
        n_allow.n(1)    = ngustpos_full(PHAA_ind);
        [n_allow.AoA(1) n_allow.Cd(1)] = find_n_conditions(ngustpos_full(PHAA_ind),...
                            v_plot(PHAA_ind),rho,S,naca2415,wgt_max);    
    end
    n_allow.V(1)    = v_plot(PHAA_ind);

    % PLAA LOAD FACTOR
    if(ngustpos_full(PLAA_ind) < npos_full(PLAA_ind))
        n_allow.name(2) = cellstr('PLAA (Maneuver)');
        n_allow.n(2)    = npos_full(PLAA_ind);
        [n_allow.AoA(2) n_allow.Cd(2)] = find_n_conditions(npos_full(PLAA_ind),...
                            v_plot(PLAA_ind),rho,S,naca2415,wgt_max);
    else
        n_allow.name(2) = cellstr('PLAA (Gust)');
        n_allow.n(2)    = ngustpos_full(PLAA_ind);
        [n_allow.AoA(2) n_allow.Cd(2)] = find_n_conditions(ngustpos_full(PLAA_ind),...
                            v_plot(PLAA_ind),rho,S,naca2415,wgt_max);        
    end
    n_allow.V(2)    = v_plot(PLAA_ind);

    % NHAA LOAD FACTOR
    if(ngustneg_full(NHAA_ind) > nneg_full(NLAA_ind))
        n_allow.name(3) = cellstr('NHAA (Maneuver)');
        n_allow.n(3)    = nneg_full(NHAA_ind);
        [n_allow.AoA(3) n_allow.Cd(3)] = find_n_conditions(nneg_full(NHAA_ind),...
                            v_plot(NHAA_ind),rho,S,naca2415,wgt_max);    
    else
        n_allow.name(3) = cellstr('NHAA (Gust)');
        n_allow.n(3)    = ngustneg_full(NHAA_ind);
        [n_allow.AoA(3) n_allow.Cd(3)] = find_n_conditions(ngustneg_full(NHAA_ind),...
                            v_plot(NHAA_ind),rho,S,naca2415,wgt_max);        
    end
    n_allow.V(3)    = v_plot(NHAA_ind);

    % NLAA LOAD FACTOR
    if(ngustneg_full(NLAA_ind) > nneg_full(NLAA_ind))
        n_allow.name(4) = cellstr('NLAA (Maneuver)');
        n_allow.n(4)    = nneg_full(NLAA_ind);
        [n_allow.AoA(4) n_allow.Cd(4)] = find_n_conditions(nneg_full(NLAA_ind),...
                            v_plot(NLAA_ind),rho,S,naca2415,wgt_max);    
    else
        n_allow.name(4) = cellstr('NLAA (Gust)');
        n_allow.n(4)    = ngustneg_full(NLAA_ind);
        [n_allow.AoA(4) n_allow.Cd(4)] = find_n_conditions(ngustneg_full(NLAA_ind),...
                            v_plot(NLAA_ind),rho,S,naca2415,wgt_max);    
    end
    n_allow.V(4)    = v_plot(NLAA_ind);

    % DIVE VELOCITY (NEGATIVE)
    if(ngustneg_full(end) > nneg_full(end))
        n_allow.name(5) = cellstr('V_{dive} Negative Load (Maneuver)');
        n_allow.n(5)    = nneg_full(end);
        [n_allow.AoA(5) n_allow.Cd(5)] = find_n_conditions(nneg_full(end),...
                    v_plot(end),rho,S,naca2415,wgt_max);    
    else
        n_allow.name(5) = cellstr('V_{dive} Negative Load (Gust)');
        n_allow.n(5)    = ngustneg_full(end);
        [n_allow.AoA(5) n_allow.Cd(5)] = find_n_conditions(ngustneg_full(end),...
                    v_plot(end),rho,S,naca2415,wgt_max);    
    end
    n_allow.V(5)    = v_plot(end);
    
    % GUST LOAD @ CRUISE
    if((ngustpos_full(NLAA_ind) > npos_full(NLAA_ind)) && ...
            (NLAA_ind > PHAA_ind))
        n_allow.name(6) = cellstr('V_{cruise} Positive Load (Gust)');
        n_allow.n(6)    = ngustpos_full(NLAA_ind);
        [n_allow.AoA(6) n_allow.Cd(6)] = find_n_conditions(ngustpos_full(NLAA_ind),...
            v_plot(NLAA_ind),rho,S,naca2415,wgt_max);    
    else
        n_allow.name(6) = cellstr('No Gust Load @ V_{cruise}');
        n_allow.n(6)    = NaN;
        n_allow.AoA(6)  = NaN;
        n_allow.Cd(6)   = NaN;
    end
    n_allow.V(6)    = v_plot(NLAA_ind);

end