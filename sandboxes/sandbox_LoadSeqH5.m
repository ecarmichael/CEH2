%% sandbox_LoadSeqH5



% cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\output_REM\seqNMF')


f_list = dir('*.h5'); 


% for iF = 1:length(f_list)
    
   this_h = MS_h5_to_stuct('seqReplayResults_LTD1_540_wake_wake.h5')
    
    
    
    
    
    
% end

%% collect the PyCaan Seq output
m_list = {'535', '537', '540'}; 
d_list = {'LTD1', 'LTD5'}; 

j_seq = []; 

cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\output_REM\seqNMF')

f_list = dir('seqReplayResults*_wake_*'); 

for ii = length(f_list):-1:1
    
    this_seq = MS_h5_to_stuct(f_list(ii).name); 
    
    m_idx = find(contains(m_list, this_seq.mouse)); 
    d_idx = find(contains(d_list, this_seq.condition));

    j_seq.(this_seq.state_ref).(this_seq.state_pred).num(d_idx, m_idx,1) = this_seq.S1_numSeqs;
    j_seq.(this_seq.state_ref).(this_seq.state_pred).num(d_idx, m_idx,2) = this_seq.S2_numSeqs;
    
        j_seq.(this_seq.state_ref).(this_seq.state_pred).rate(d_idx, m_idx,1) = this_seq.S1_numSeqs/length(this_seq.H)/30;
    j_seq.(this_seq.state_ref).(this_seq.state_pred).rate(d_idx, m_idx,2) = this_seq.S2_numSeqs/length(this_seq.H)/30;
    
    
    fprintf( '%s %s   N %0.0f   R: %0.2f  dur: %0.2f\n', this_seq.state_ref, this_seq.state_pred, this_seq.S1_numSeqs, this_seq.S1_numSeqs/length(this_seq.H)/30, length(this_seq.H)/30)
end





