function pos = MS_DLC_cat(cfg_in, DLC_dir)
%% MS_DLC_cat: load and merge multiple DLC outputs. 
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-09-23   initial version 
%
%
%
%% initialize

cfg_def = [];
cfg_def.conv_fac = [1 1];
cfg_def.model_in = []; 
cfg_def.plot_flag = []; 

cfg = ProcessConfig(cfg_def, cfg_in); 

%%  Get the position data if it is there


    
    if length(DLC_dir) <2
    
        [this_pos, ~] = MS_DLC2TSD(DLC_dir, [], cfg.conv_fac);
%         if pos.tvec(1) <0
%             pos.tvec = pos.tvec + abs(pos.tvec(1));
%         elseif pos.tvec(1) >0
%             pos.tvec  = pos.tvec - abs(pos.tvec(1));
%         end
    
    else
        this_pos = []; 
        all_tvec = []; all_data = [];
        for ii = 1:length(DLC_dir)
            [this_pos, ~] = MS_DLC2TSD(DLC_dir{ii}, cfg.model_in, cfg.conv_fac, cfg.plot_flag);

            
            if ii == 1
                all_data = this_pos.data; 
                all_tvec =  this_pos.tvec; 
            else
                all_tvec =  [all_tvec;  this_pos.tvec];
                all_data = [all_data, this_pos.data]; 

            end
%             all_tvec = [all_tvec, this_pos{ii}.tvec']; 
%             all_data = [all_data, this_pos{ii}.data];
%             this_pos{1} = tsd(all_tvec', all_data, 'label',this_pos{1}.label);
%          this_pos{1}.units = this_pos{1}.units; 
            fprintf('<strong>%s</strong>: Recoding #%.0f: duration = %.2fmin\n', mfilename,ii,  (this_pos.tvec(end) - this_pos.tvec(1))/60)

        end
        
        % merge the pos files
%         pos = tsd(all_tvec, all_data, 'label',this_pos{1}.label);
%         pos.units = this_pos{1}.units; 
        
    end


pos = this_pos; 
pos.tvec = all_tvec; 
pos.data = all_data; 