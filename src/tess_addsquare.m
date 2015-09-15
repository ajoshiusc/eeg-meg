function [TessMat, errMsg] = tess_addsquare(TessFile, SquareFile, AtlasSquareFile)
% TESS_ADD: Add a BrainSuite registered square to an existing surface.
%
% USAGE:  TessMat = tess_addSquare(TessFile, SquareFile=select)
%         TessMat = tess_addSquare(TessMat,  SquareFile=select)

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2013

% Initialize returned variables
TessMat = [];
errMsg = [];

% Ask for Square file
if (nargin < 3) || isempty(SquareFile) || isempty(AtlasSquareFile)
    % Get last used directories and formats
    LastUsedDirs = bst_get('LastUsedDirs');
    % Get Surface files
    SquareFile = java_getfile( 'open', ...
       'Import atlas surfaces...', ...      % Window title
       LastUsedDirs.ImportAnat, ...   % Default directory
       'single', 'files', ...         % Selection mode
       {{'.dfs'}, 'Registered SVReg Square (*.dfs)', 'BrainSuite'}, 'BrainSuite');
    % If no file was selected: exit
    if isempty(AtlasSquareFile)
        return
    end
    % Save default import directory
    LastUsedDirs.ImportAnat = bst_fileparts(AtlasSquareFile);
    bst_set('LastUsedDirs', LastUsedDirs);
end

% Progress bar
isProgressBar = ~bst_progress('isVisible');
if isProgressBar
    bst_progress('start', 'Load registration', 'Loading BrainSuite Square...');
end

% Get the subject MRI
[sSubject, iSubject] = bst_get('SurfaceFile', TessFile);
if isempty(sSubject.Anatomy) || isempty(sSubject.Anatomy(1).FileName)
    errMsg = 'Subject does not have a registered MRI.';
    return;
end
sMri = bst_memory('LoadMri', iSubject);
    
% If destination surface is already loaded
if isstruct(TessFile)
    TessMat = TessFile;
    TessFile = [];
% Else: load target surface file
else
    TessMat = in_tess_bst(TessFile);
end

% Load the Square surface: DO NOT CONVERT TO SCS!!!!
%SquareMat = in_tess(SquareFile, 'FS', sMri);
% Load the surface, keep in the original coordinate system
surfSqr = readdfs(SquareFile);
atlasSqr= readdfs(AtlasSquareFile);

SquareVertices=[surfSqr.u',surfSqr.v'];
AtlasSquareVertices=[atlasSqr.u',atlasSqr.v'];
% Check that the number of vertices match
if (length(SquareVertices) ~= length(TessMat.Vertices))
    errMsg = sprintf('The number of vertices in the surface (%d) and the Square (%d) do not match.', length(TessMat.Vertices), length(SquareMat.Vertices));
    TessMat = [];
    return;
end
% Add the Square vertex information to the surface matrix
TessMat.Reg.Square.Vertices = SquareVertices;
TessMat.Reg.AtlasSquare.Vertices = AtlasSquareVertices;
% Save modifications to input file
bst_save(file_fullpath(TessFile), TessMat, 'v7');

% Close progress bar
if isProgressBar
    bst_progress('stop');
end





