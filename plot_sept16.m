%%
%Parameters
channel = 'GA_64';
ch_idx = 48;
event_type = {'A','V','BP1','BP2','BP3','BP4','BP5'};
%%
%Event type 1, Auditory word
A_end = find(.75<hgp_trials_comb.timeavg,1);
A_dat = hgp_trials_comb.avgdat{1}.avg(48,1:A_end-1);

%%
%Event type 2, Visual word
V_end = find(1.25<hgp_trials_comb.timeavg,1);
V_dat = hgp_trials_comb.avgdat{1}.avg(48,A_end:V_end-1);

%%
%Event type 3, Button press 1
BP1_avgrt = 2.4386;
BP1_start = find(BP1_avgrt-.2<hgp_trials_comb.timeavg,1);
BP1_end = find(BP1_avgrt+.1<hgp_trials_comb.timeavg,1);
BP1_dat = hgp_trials_comb.avgdat{1}.avg(48,BP1_start:BP1_end);

%%
%Event type 4, Button press 2

BP2_avgrt = 3.9279;
BP2_start = find(BP2_avgrt-.2<hgp_trials_comb.timeavg,1);
BP2_end = find(BP2_avgrt+.1<hgp_trials_comb.timeavg,1);
BP2_dat = hgp_trials_comb.avgdat{1}.avg(48,BP2_start:BP2_end);

%%
%Event type 5, Button press 3

BP3_avgrt = 5.414;
BP3_start = find(BP3_avgrt-.2<hgp_trials_comb.timeavg,1);
BP3_end = find(BP3_avgrt+.1<hgp_trials_comb.timeavg,1);
BP3_dat = hgp_trials_comb.avgdat{1}.avg(48,BP3_start:BP3_end);

%%
%Event type 6, Button press 4

BP4_avgrt = 6.8618;
BP4_start = find(BP4_avgrt-.2<hgp_trials_comb.timeavg,1);
BP4_end = find(BP4_avgrt+.1<hgp_trials_comb.timeavg,1);
BP4_dat = hgp_trials_comb.avgdat{1}.avg(48,BP4_start:BP4_end);

%%
%Event Type 7, Button press 5

BP5_avgrt = 8.3283;
BP5_start = find(BP5_avgrt-.2<hgp_trials_comb.timeavg,1);
BP5_end = find(BP5_avgrt+.1<hgp_trials_comb.timeavg,1);
BP5_dat = hgp_trials_comb.avgdat{1}.avg(48,BP5_start:BP5_end);


%%
figure; hold on;
plot(hgp_trials_comb.timeavg(0:A_end),A_dat);hold on;
set(gca,'XTick',[0 A_end/2 A_end],'XTickLabel',{'0' 'halfaud' 'fullaud'});
%%
% figure; hold on;
plot(V_dat);hold on;
set(gca,'XTick',[A_end V_end/2 V_end],'XTickLabel',{'A_end' 'halfv' 'fullv'});

% plot(V_dat);hold on;
% plot(BP1_dat);hold on;
%set(gca,'XTick',[0 .5 .75],'XTickLabel',{'0' '500' '750'});
%%
%Create plot
for i = 1 : length(event_type) % per event type
    