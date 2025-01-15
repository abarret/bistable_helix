function [avg_vals, max_vals, min_vals] = getAvgPitchAndRadius(prs, t)
t_vals = prs(1,:);
t_idxs = t_vals > t;
avg_vals = [mean(prs(2,t_idxs)), mean(prs(3,t_idxs))];
max_vals = [max(prs(2,t_idxs)), max(prs(3,t_idxs))];
min_vals = [min(prs(2,t_idxs)), min(prs(3,t_idxs))];
end