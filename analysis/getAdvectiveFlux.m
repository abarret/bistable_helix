function fl = getAdvectiveFlux(base_path)
%% Get the advective flux file
file = base_path + "/advective_flux.txt";
%% Read file
% File contains time and flux
fileID = fopen(file, 'r');
fl = fscanf(fileID, '%f %f', [2, inf]);
fclose(fileID);
end