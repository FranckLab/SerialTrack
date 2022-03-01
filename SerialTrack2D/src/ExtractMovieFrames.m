% Demo to extract frames and get frame means from a movie and save individual frames to separate image files.
% Then rebuilds a new movie by recalling the saved images from disk.
% Also computes the mean gray value of the color channels
% And detects the difference between a frame and the previous frame.
% Illustrates the use of the VideoReader and VideoWriter classes.
% A Mathworks demo (different than mine) is located here 
% http://www.mathworks.com/help/matlab/examples/convert-between-image-sequences-and-video.html

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 22;

% Open the rhino.avi demo movie that ships with MATLAB.
% First get the folder that it lives in.
videoname = '1242248s2.avi';
folder = fileparts(which(videoname)); % Determine where demo folder is (works with all versions).
% Pick one of the two demo movies shipped with the Image Processing Toolbox.
% Comment out the other one.
movieFullFileName = fullfile(folder, videoname);
% movieFullFileName = fullfile(folder, 'traffic.avi');
% Check to see that it exists.
if ~exist(movieFullFileName, 'file')
	strErrorMessage = sprintf('File not found:\n%s\nYou can choose a new one, or cancel', movieFullFileName);
	response = questdlg(strErrorMessage, 'File not found', 'OK - choose a new movie.', 'Cancel', 'OK - choose a new movie.');
	if strcmpi(response, 'OK - choose a new movie.')
		[baseFileName, folderName, FilterIndex] = uigetfile('*.avi');
		if ~isequal(baseFileName, 0)
			movieFullFileName = fullfile(folderName, baseFileName);
		else
			return;
		end
	else
		return;
	end
end


%% 
%     numberOfFrames = 526;
%     for frame = 1: numberOfFrames
%         thisFrame =  double(rgb2gray(read(videoObject, frame)));
%         if frame == 1
%             BackgroundFrame = thisFrame;
%         else
%             BackgroundFrame = BackgroundFrame + thisFrame;
%         end  
%     end
%     BackgroundFrame =  double(BackgroundFrame)/numberOfFrames;
%     figure, subplot(2,2,1); imshow(uint8(BackgroundFrame));
%     
%     frame = 1; firstFrame =  double(rgb2gray(read(videoObject, frame)));
%     subplot(2,2,2); imshow(uint8(firstFrame));
%     frame = 2; secondFrame =  double(rgb2gray(read(videoObject, frame)));
%     subplot(2,2,3); imshow(uint8(secondFrame));
%     subplot(2,2,4); imshow(uint8(firstFrame-BackgroundFrame)); caxis auto;
    
    
%%
try
	videoObject = VideoReader(movieFullFileName)
	% Determine how many frames there are.
	numberOfFrames = videoObject.NumberOfFrames;
	vidHeight = videoObject.Height;
	vidWidth = videoObject.Width;
	
	numberOfFramesWritten = 0;
	% Prepare a figure to show the images in the upper half of the screen.
	figure;
	% 	screenSize = get(0, 'ScreenSize');
	% Enlarge figure to full screen.
	set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
	
	% Ask user if they want to write the individual frames out to disk.
	promptMessage = sprintf('Do you want to save the individual frames out to individual disk files?');
	button = questdlg(promptMessage, 'Save individual frames?', 'Yes', 'No', 'Yes');
	if strcmp(button, 'Yes')
		writeToDisk = true;
		
		% Extract out the various parts of the filename.
		[folder, baseFileName, extentions] = fileparts(movieFullFileName);
		% Make up a special new output subfolder for all the separate
		% movie frames that we're going to extract and save to disk.
		% (Don't worry - windows can handle forward slashes in the folder name.)
		folder = pwd;   % Make it a subfolder of the folder where this m-file lives.
		outputFolder = sprintf('%s/Movie Frames from %s', folder, baseFileName);
		% Create the folder if it doesn't exist already.
		if ~exist(outputFolder, 'dir')
			mkdir(outputFolder);
		end
	else
		writeToDisk = false;
	end
	
	% Loop through the movie, writing all frames out.
	% Each frame will be in a separate file with unique name.
	meanGrayLevels = zeros(numberOfFrames, 1);
	meanRedLevels = zeros(numberOfFrames, 1);
	meanGreenLevels = zeros(numberOfFrames, 1);
	meanBlueLevels = zeros(numberOfFrames, 1);
    thresholdLevels = zeros(numberOfFrames, 1);
    
    
    
	for frame = 1 : numberOfFrames
		% Extract the frame from the movie structure.
		thisFrame = rgb2gray(read(videoObject, frame));
		
		% Display it
		hImage = subplot(2, 2, 1);
		imshow(thisFrame);
		caption = sprintf('Frame %4d of %d.', frame, numberOfFrames);
		title(caption, 'FontSize', fontSize);
		drawnow; % Force it to refresh the window.
		
		% Write the image array to the output file, if requested.
		if writeToDisk
			% Construct an output image file name.
			outputBaseFileName = sprintf('Frame %4.4d.tiff', frame);
			outputFullFileName = fullfile(outputFolder, outputBaseFileName);
			
			% Stamp the name and frame number onto the image.
			% At this point it's just going into the overlay,
			% not actually getting written into the pixel values.
			% JY!!! % text(5, 15, outputBaseFileName, 'FontSize', 20);
			
			% Extract the image with the text "burned into" it.
			frameWithText = getframe(gca);
			% frameWithText.cdata is the image with the text
			% actually written into the pixel values.
			% Write it out to disk.
			imwrite((thisFrame), outputFullFileName, 'tiff');
		end
		
		% Calculate the mean gray level.
		grayImage = (thisFrame);
		meanGrayLevels(frame) = mean(grayImage(:));
		
		% % Calculate the mean R, G, and B levels.
		%meanRedLevels(frame) = mean(mean(thisFrame(:, :, 1)));
		%meanGreenLevels(frame) = mean(mean(thisFrame(:, :, 2)));
		%meanBlueLevels(frame) = mean(mean(thisFrame(:, :, 3)));
		
		% Plot the mean gray levels.
		hPlot = subplot(2, 2, 2);
		%hold off;
		%plot(meanGrayLevels, 'k-', 'LineWidth', 3);
		%hold on;
		%plot(meanRedLevels, 'r-', 'LineWidth', 2);
		%plot(meanGreenLevels, 'g-', 'LineWidth', 2);
		%plot(meanBlueLevels, 'b-', 'LineWidth', 2);
		%grid on;
		
		% Put title back because plot() erases the existing title.
		title('Mean Gray Levels', 'FontSize', fontSize);
		if frame == 1
			xlabel('Frame Number');
			ylabel('Gray Level');
			% Get size data later for preallocation if we read
			% the movie back in from disk.
			[rows, columns, numberOfColorChannels] = size(thisFrame);
		end
		
		% Update user with the progress.  Display in the command window.
		if writeToDisk
			progressIndication = sprintf('Wrote frame %4d of %d.', frame, numberOfFrames);
		else
			progressIndication = sprintf('Processed frame %4d of %d.', frame, numberOfFrames);
		end
		disp(progressIndication);
		% Increment frame count (should eventually = numberOfFrames
		% unless an error happens).
		numberOfFramesWritten = numberOfFramesWritten + 1;
		
% 		% Now let's do the differencing
% 		alpha = 1;
% 		if frame == 1
% 			%Background = thisFrame;
%             Background = (1-alpha)*double(thisFrame) + alpha*double(BackgroundFrame);
% 		else
% 			% Change background slightly at each frame
% 			% 			Background(t+1)=(1-alpha)*I+alpha*Background
% 			%Background = (1-alpha)* thisFrame + alpha * Background;
%             Background = (1-alpha)*double(thisFrame) + alpha*double(BackgroundFrame);
% 		end
% 		% Display the changing/adapting background.
% 		subplot(2, 2, 3);
% 		imshow(Background);
% 		title('Adaptive Background', 'FontSize', fontSize);
% 		% Calculate a difference between this frame and the background.
% 		differenceImage = thisFrame - uint8(Background);
% 		% Threshold with Otsu method.
% 		grayImage =  (differenceImage); % Convert to gray level
% 		thresholdLevel = 0.1; % graythresh(grayImage); % Get threshold.
%         thresholdLevels(frame) = thresholdLevel;
%         subplot(2, 2, 2); 
% 		plot(thresholdLevels, 'k-', 'LineWidth', 2);
% 		binaryImage = im2bw( grayImage, thresholdLevel); % Do the binarization
% 		% Plot the binary image.
% 		subplot(2, 2, 4);
% 		imshow(binaryImage);
% 		title('Binarized Difference Image', 'FontSize', fontSize);
	end
	
	% Alert user that we're done.
	if writeToDisk
		finishedMessage = sprintf('Done!  It wrote %d frames to folder\n"%s"', numberOfFramesWritten, outputFolder);
	else
		finishedMessage = sprintf('Done!  It processed %d frames of\n"%s"', numberOfFramesWritten, movieFullFileName);
	end
	disp(finishedMessage); % Write to command window.
	uiwait(msgbox(finishedMessage)); % Also pop up a message box.
	
	% Exit if they didn't write any individual frames out to disk.
	if ~writeToDisk
		return;
	end
	
	% Ask user if they want to read the individual frames from the disk,
	% that they just wrote out, back into a movie and display it.
	promptMessage = sprintf('Do you want to recall the individual frames\nback from disk into a movie?\n(This will take several seconds.)');
	button = questdlg(promptMessage, 'Recall Movie?', 'Yes', 'No', 'Yes');
	if strcmp(button, 'No')
		return;
	end

	% Create a VideoWriter object to write the video out to a new, different file.
	writerObj = VideoWriter('NewRhinos.avi');
	open(writerObj);
	
	% Read the frames back in from disk, and convert them to a movie.
	% Preallocate recalledMovie, which will be an array of structures.
	% First get a cell array with all the frames.
	allTheFrames = cell(numberOfFrames,1);
	allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
	% Next get a cell array with all the colormaps.
	allTheColorMaps = cell(numberOfFrames,1);
	allTheColorMaps(:) = {zeros(256, 3)};
	% Now combine these to make the array of structures.
	recalledMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps)
	for frame = 1 : numberOfFrames
		% Construct an output image file name.
		outputBaseFileName = sprintf('Frame %4.4d.png', frame);
		outputFullFileName = fullfile(outputFolder, outputBaseFileName);
		% Read the image in from disk.
		thisFrame = imread(outputFullFileName);
		% Convert the image into a "movie frame" structure.
		recalledMovie(frame) = im2frame(thisFrame);
		% Write this frame out to a new video file.
		writeVideo(writerObj, thisFrame);
	end
	close(writerObj);
	% Get rid of old image and plot.
	delete(hImage);
	delete(hPlot);
	% Create new axes for our movie.
	subplot(1, 3, 2);
	axis off;  % Turn off axes numbers.
	title('Movie recalled from disk', 'FontSize', fontSize);
	% Play the movie in the axes.
	movie(recalledMovie);
	% Note: if you want to display graphics or text in the overlay
	% as the movie plays back then you need to do it like I did at first
	% (at the top of this file where you extract and imshow a frame at a time.)
	msgbox('Done with this demo!');
	
catch ME
	% Some error happened if you get here.
	strErrorMessage = sprintf('Error extracting movie frames from:\n\n%s\n\nError: %s\n\n)', movieFullFileName, ME.message);
	uiwait(msgbox(strErrorMessage));
end

