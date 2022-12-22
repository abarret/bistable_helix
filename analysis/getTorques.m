function torques = getTorques(base_path)
%% Get list of files
sys_files = dir(base_path);
files = cell(1, length(sys_files)-2);
% Skip first two since those are . and ..
for i = 1:length(sys_files)-2
    files{i} = sys_files(i+2).name;
end
%% Read files
% File contains three torque components.
% Time is given by numbers before .out
dt = 1.0e-7;
torques = zeros(4,length(files));
for i = 1:length(files)
    % Determine time
    time = dt * str2double(regexp(files{i}, '\d+', 'match', 'once'));
    fileID = fopen(base_path + "/" + files{i}, 'r');
    % Should only be one line. Read it.
    data = fscanf(fileID, '%f %f %f');
    fclose(fileID);
    torques(1,i) = time;
    torques(2:end,i) = data;
end

%% Sort according to time
times = torques(1,:);
[~,idxs] = sort(times);
torques = torques(:,idxs);

end