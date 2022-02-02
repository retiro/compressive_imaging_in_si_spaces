%RESULTSSETUP Setup file for the Matrix-free IPM solver.
%
%  This file installs the appropriate software packages for the Matrix-free IPM
%
%  Actions
%  -------
%
%    Adds to path the files of the current folder:
%
% Copyright 2013, Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich

% Get root location of resultsSetup
root = fileparts(which(mfilename)); 

% ----------------------------------------------------------------------
% Add the appropriate subdirs to path.
% ----------------------------------------------------------------------

addpath(genpath(root));
fprintf('%s \n',root);

fprintf(['\nThe above directories and their subdirectories have been'  ,...
         ' temporarily added to your path.\n']);

% ---------------------------------------------------------------------
% Store the new path permanently.
% ---------------------------------------------------------------------
while (true)
    reply = input('Would you like to save the new path (optional)? Y/N [N]: ', 's');
    if (strcmpi(reply,'Y') || strcmpi(reply,'N'))
        break;   
    end
end

if (strcmpi(reply,'Y'))
    is_saved = savepath;
    if (~is_saved)
        disp('The new path is saved.');
    else
        startfold = userpath;
        fprintf('\n%s \n \n','The new path could not be saved.');
        fprintf('Please copy the following addpath commands in your\n');
        fprintf('startup.m file.\n\n');
        text = ['addpath(''',root,''')'];
        fprintf('%s \n\n',text);
        fprintf('The startup.m file is located in MATLAB startup folder but it is not by default created.\n');
        fprintf('It has to be created by the user. For more information see:\n');
        fprintf('http://www.mathworks.co.uk/help/matlab/ref/startup.html\n\n');
        fprintf(['Your MATLAB startup folder is located at:\n\n',startfold,'\n\n']);
        fprintf('For information about MATLAB startup folder see:\n');
        fprintf('http://www.mathworks.co.uk/help/matlab/matlab_env/matlab-startup-folder.html\n');
    end
else
    disp('The new path was not saved.');
end

clear;
