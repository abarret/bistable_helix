function fl = getAdvectiveFlux(base_path)
%% Get the KE file
file = base_path + "/ke.txt";
%% Read file
% File contains time and KE
fileID = fopen(file, 'r');
fl = fscanf(fileID, '%f %f', [2, inf]);
fclose(fileID);
end