function dSub_get_behav(evts)
%% dSub_get_behav:  collects the trial data for the W maze experiments.
%
%
%
%     Inputs
%        - Meta: [struct]   Optional meta data generated with
%        MS_Write_meta_dSub
%
%
%
%
%% Get Meta data if it exists

% if nargin < 1
%     try
%         MS_Load_meta
%     catch
%         try
%             MS_Write_meta_dSub;
%         catch
%             error('Can''t load or generate Meta data.')
%         end
%     end
% end
%
%
%
% %% Get the events



%% Split out entries and exits

b_ent_idx = contains(evt.label, 'Zoned Video: Box Entered');
b_ext_idx = contains(evt.label, 'Zoned Video: Box Exited');

c_ent_idx = contains(evt.label, 'Zoned Video: Center Entered');
c_ext_idx = contains(evt.label, 'Zoned Video: Center Exited');

r_ent_idx = contains(evt.label, 'Zoned Video: R_reward Entered');
r_ext_idx = contains(evt.label, 'Zoned Video: R_reward Exited');

l_ent_idx = contains(evt.label, 'Zoned Video: L_reward Entered');
l_ext_idx = contains(evt.label, 'Zoned Video: L_reward Exited');

b_idx = contains(evt.label, 'b');
s_idx = contains(evt.label, 's');

r_idx = contains(evt.label, 'r');
l_idx = contains(evt.label, 'l');



b_ent_evt = evt.t{b_ent_idx};
b_ext_evt = evt.t{b_ext_idx};

c_ent_evt = evt.t{c_ent_idx};
c_ext_evt = evt.t{c_ext_idx};

r_ent_evt = evt.t{r_ent_idx};
r_ext_evt = evt.t{r_ext_idx};

l_ent_evt = evt.t{l_ent_idx};
l_ext_evt = evt.t{l_ext_idx};

b_evt = evt.t{b_idx};
s_evt = evt.t{s_idx};

r_evt = evt.t{r_idx};
l_evt = evt.t{l_idx};


%% get the specific times for each event type.

% get the box to center 'trial start'
b2c = [];
for iC = 1:length(c_ent_evt)
    if ~isempty(find(diff(b_ext_evt < c_ent_evt(iC)) == -1))
        b2c(iC) = find(diff(b_ext_evt < c_ent_evt(iC)) == -1);
    end
end
if ~isempty(b2c); b2c = unique(b2c); b2c_t = evt.t{c_ent_idx}(b2c); end

% get the Reward to box entries
    r2b = []; l2b = [];
for iB = 1:length(b_ent_evt)
    if ~isempty(find(diff(r_ent_evt < b_ent_evt(iB)) == -1))
        r2b(iB) = find(diff(r_ent_evt < b_ent_evt(iB)) == -1);
    end
    if ~isempty(find(diff(r_ent_evt < b_ent_evt(iB)) == -1))
        l2b(iB) = find(diff(l_ent_evt < b_ent_evt(iB)) == -1);
    end
end
if ~isempty(r2b); r2b = unique(r2b); r2b_t = evt.t{b_ent_idx}(r2b); end
if ~isempty(l2b); l2b = unique(l2b); l2b_t = evt.t{b_ent_idx}(l2b); end



% get the Reward t box
    r2b = []; l2b = [];
for iB = 1:length(b_ent_evt)
    if ~isempty(find(diff(r_ent_evt < b_ent_evt(iB)) == -1))
        r2b(iB) = find(diff(r_ent_evt < b_ent_evt(iB)) == -1);
    end
    if ~isempty(find(diff(r_ent_evt < b_ent_evt(iB)) == -1))
        l2b(iB) = find(diff(l_ent_evt < b_ent_evt(iB)) == -1);
    end
end
if ~isempty(r2b); r2b = unique(r2b); r2b_t = evt.t{b_ent_idx}(r2b); end
if ~isempty(l2b); l2b = unique(l2b); l2b_t = evt.t{b_ent_idx}(l2b); end




