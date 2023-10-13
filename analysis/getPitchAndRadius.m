function pr = getPitchAndRadius(base_path, idxs)
%% Get list of files
sys_files = dir(base_path);
files = cell(1, length(sys_files)-2);
% Skip first two since those are . and ..
for i = 1:length(sys_files)-2
    files{i} = sys_files(i+2).name;
end
%% Read files
% File contains index, radius, and pitch.
% Time is given by numbers before .out
dt = 1.0e-7;
pr = zeros(3, length(files),length(idxs));
for i = 1:length(files)
    % Determine time
    %time = dt * str2double(regexp(files{i}, '\d+', 'match', 'once'));
    fileID = fopen(base_path + "/" + files{i}, 'r');
    % Should be many lines consisting of three elements. Read it.
    time = str2double(fgetl(fileID));
    data = fscanf(fileID, '%f %f %f', [3, inf]);
    % Grab 
    fclose(fileID);
    
    for j = 1:length(idxs)
        pr(1,i,j) = time;
        pr(2:3,i,j) = data(2:3,idxs(j));
    end
end

%% Sort according to time
times = pr(1,:,1);
[~,idxs] = sort(times);
pr = pr(:,idxs,:);
end