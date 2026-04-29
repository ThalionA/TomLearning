%% Aggregate_Sweep_Results_to_PPT.m
import mlreportgen.ppt.*

% Define the target directories
folders = {
    '/Users/theoamvr/Desktop/Experiments/TomLearning/HC_V1_data/Sweep_Results_2026_03_02',
    '/Users/theoamvr/Desktop/Experiments/TomLearning/HC_V1_data/Sweep_Results_2026_02_27'
};

% Define the output PowerPoint path
pptPath = '/Users/theoamvr/Desktop/Experiments/TomLearning/HC_V1_data/Epoch_Precession_Summary.pptx';

% Initialize Presentation
ppt = Presentation(pptPath);
open(ppt);

% Add a title slide
titleSlide = add(ppt, 'Title Slide');
replace(titleSlide, 'Title', 'Epoch Precession Results');
replace(titleSlide, 'Subtitle', 'Hyperparameter Sweep Aggregation');

count = 0;

% Iterate through folders and extract figures
for f = 1:length(folders)
    currentFolder = folders{f};
    
    if ~exist(currentFolder, 'dir')
        fprintf('Directory not found, skipping: %s\n', currentFolder);
        continue;
    end
    
    % Find both SVG and PNG (fallback) files matching the epoch precession prefix
    files = [dir(fullfile(currentFolder, 'Epoch_Precession_*.svg')); ...
             dir(fullfile(currentFolder, 'Epoch_Precession_*.png'))];
         
    for i = 1:length(files)
        fileName = files(i).name;
        filePath = fullfile(currentFolder, fileName);
        
        % Create a new slide with Title and Content layout
        slide = add(ppt, 'Title and Content');
        
        % Format the slide title (remove prefix and extension)
        slideTitle = strrep(fileName, 'Epoch_Precession_', '');
        slideTitle = regexprep(slideTitle, '\.(svg|png)$', '');
        
        replace(slide, 'Title', slideTitle);
        
        % Insert the figure
        img = Picture(filePath);
        replace(slide, 'Content', img);
        
        count = count + 1;
    end
end

% Close and save the presentation
close(ppt);

fprintf('Successfully added %d figures to the presentation.\n', count);
fprintf('Presentation saved to: %s\n', pptPath);

% Automatically open the generated file
if ismac
    system(['open "', pptPath, '"']);
end