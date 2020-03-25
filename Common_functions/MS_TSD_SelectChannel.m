function tsd_out = MS_TSD_SelectChannel(tsd_in,channel_label)
% function tsd_out = TSD_SelectChannel(tsd_in,channel_label)
%
% returns tsd containing only channels with labels in channel_label
%
% MvdM 2015-11-03 initial version
% EC 2020-03-20 added removal of cfg.hdr cells.  

if ~iscell(channel_label)
   keep_idx = strmatch(channel_label,tsd_in.label,'exact');
else
   keep_idx = [];
   for iCh = 1:length(channel_label)
       
      keep_idx = cat(1,keep_idx,strmatch(channel_label{iCh},tsd_in.label,'exact'));
       
   end
end

if isempty(keep_idx)
    error('No channels selected.');
end

tsd_out = tsd_in;
tsd_out.data = tsd_out.data(keep_idx,:);
tsd_out.label = tsd_out.label(keep_idx);
tsd_out.cfg.hdr = tsd_out.cfg.hdr(keep_idx); 