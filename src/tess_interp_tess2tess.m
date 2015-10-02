function [Wmat, sSrcSubj, sDestSubj, srcSurfMat, destSurfMat, isStopWarped] = tess_interp_tess2tess( srcSurfFile, destSurfFile, isInteractive, isStopWarped )
% TESS_INTERP_TESS2TESS: Compute an inteprolation matrix between two cortex surfaces.

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
% Authors: Francois Tadel, Anand Joshi in 2015

% Parse inputs
if (nargin < 4) || isempty(isStopWarped)
    isStopWarped = [];
end
if (nargin < 3) || isempty(isInteractive)
    isInteractive = 1;
end

% ===== GET SURFACES =====
% Load source surface file
srcSurfMat  = in_tess_bst(srcSurfFile);
% Load destination surface file
destSurfMat = in_tess_bst(destSurfFile);
% Get source and destination subjects
sSrcSubj  = bst_get('SurfaceFile', srcSurfFile);
sDestSubj = bst_get('SurfaceFile', destSurfFile);
% Number of vertices
nSrc  = size(srcSurfMat.Vertices, 1);
nDest = size(destSurfMat.Vertices, 1);
% Source subject and destination subject are the same
isSameSubject = file_compare(sSrcSubj.FileName, sDestSubj.FileName);
% Check if source or destination are the default anatomy
isSrcDefaultSubj  = ismember(bst_fileparts(sSrcSubj.FileName),  {bst_get('DirDefaultSubject'), bst_get('NormalizedSubjectName')});
isDestDefaultSubj = ismember(bst_fileparts(sDestSubj.FileName), {bst_get('DirDefaultSubject'), bst_get('NormalizedSubjectName')});
% Signature string for the current transformation
Signature = sprintf('%s%d=>%s%d', srcSurfFile, length(srcSurfMat.Vertices), destSurfFile, length(destSurfMat.Vertices));
% Initialize interpolation matrix
Wmat = [];
isSaveInterp = 1;


% ===== RE-USE PREVIOUS INTERPOLATION =====
% Try to get an existing valid interpolation matrix
if isempty(Wmat) && isfield(srcSurfMat, 'tess2tess_interp') && all(isfield(srcSurfMat.tess2tess_interp, {'Signature', 'Wmat'})) && ...
        strcmpi(srcSurfMat.tess2tess_interp.Signature, Signature) && ~isempty(srcSurfMat.tess2tess_interp.Wmat)
    Wmat = srcSurfMat.tess2tess_interp.Wmat;
    % Do not save this interpolation matrix, it's already saved
    isSaveInterp = 0;
end

% ===== CHECK IF WARPED =====
% If projecting a warped subject back on the original brain: NOT necessary
if isempty(Wmat) && ~isempty(strfind(srcSurfFile, '_warped')) && ~isSrcDefaultSubj && isDestDefaultSubj && (nSrc == nDest)
    % Warning message
    warnMsg = ['The source files were computed on a warped anatomy, there is' 10 ...
        'no need to re-project them on the default anatomy, you can directly' 10 ...
        'calculate average or differences across subjects.'];
    % Ask user to cancel the process
    if isempty(isStopWarped)
        if isInteractive
            isStopWarped = ~java_dialog('confirm', [warnMsg 10 10 'Project sources anyways?'], 'Project sources');
            if isStopWarped
                bst_progress('stop');
                return;
            end
        elseif ~isInteractive
            isStopWarped = 0;
            bst_report('Warning', 'process_project_sources', [], warnMsg);
        end
    end
    % Interpolation matrix: Use an identity matrix
    Wmat = speye(nDest,nSrc);
    % Do not save this interpolation matrix, it's really not necessary
    isSaveInterp = 0;
end

% ===== USE FREESURFER SPHERES =====
% If the registered spheres are available in both surfaces
if isempty(Wmat) && isfield(srcSurfMat, 'Reg') && isfield(srcSurfMat.Reg, 'Sphere') && isfield(srcSurfMat.Reg.Sphere, 'Vertices') && ~isempty(srcSurfMat.Reg.Sphere.Vertices) && ...
        isfield(destSurfMat, 'Reg') && isfield(destSurfMat.Reg, 'Sphere') && isfield(destSurfMat.Reg.Sphere, 'Vertices') && ~isempty(destSurfMat.Reg.Sphere.Vertices)
    % Evaluate number of vertices to use
    nbNeighbors = 8;
    % Allocate interpolation matrix
    Wmat = spalloc(nDest, nSrc, nbNeighbors * nDest);
    % Split hemispheres
    [rHsrc, lHsrc, isConnected(1)]  = tess_hemisplit(srcSurfMat);
    [rHdest,lHdest, isConnected(2)] = tess_hemisplit(destSurfMat);
    % Get vertices
    srcVert  = double(srcSurfMat.Reg.Sphere.Vertices);
    destVert = double(destSurfMat.Reg.Sphere.Vertices);
    % If hemispheres are connected: process all at once
    if any(isConnected)
        rHsrc  = 1:nSrc;
        rHdest = 1:nDest;
        lHsrc  = [];
        lHdest = [];
    end
    % Re-interpolate using the sphere and the shepards algorithm
    WmatTmp = bst_shepards(destVert(rHdest,:), srcVert(rHsrc,:), nbNeighbors, 0);
    Wmat(rHdest,rHsrc) = WmatTmp;
    if ~isempty(lHdest)
        WmatTmp = bst_shepards(destVert(lHdest,:), srcVert(lHsrc,:), nbNeighbors, 0);
        Wmat(lHdest,lHsrc) = WmatTmp;
    end
end

% ===== USE BrainSuite Squares =====
% If the registered squares are available in both surfaces
if isempty(Wmat) && isfield(srcSurfMat, 'Reg') && isfield(srcSurfMat.Reg, 'Square') && isfield(srcSurfMat.Reg.Square, 'Vertices') && ~isempty(srcSurfMat.Reg.Square.Vertices) && ...
        isfield(destSurfMat, 'Reg') && isfield(destSurfMat.Reg, 'Square') && isfield(destSurfMat.Reg.Square, 'Vertices') && ~isempty(destSurfMat.Reg.Square.Vertices)
    % Evaluate number of vertices to use
    nbNeighbors = 8;
    % Allocate interpolation matrix
    Wmat = spalloc(nDest, nSrc, nbNeighbors * nDest);
    % Split hemispheres
    [rHsrc, lHsrc, isConnected(1)]  = tess_hemisplit(srcSurfMat);
    [rHdest,lHdest, isConnected(2)] = tess_hemisplit(destSurfMat);
    % Get vertices
    srcVert  = double(srcSurfMat.Reg.Square.Vertices);
    srcAtlasVert  = double(srcSurfMat.Reg.AtlasSquare.Vertices);
    destVert = double(destSurfMat.Reg.Square.Vertices);
    destAtlasVert = double(destSurfMat.Reg.AtlasSquare.Vertices);
    
    AtlLInd=find(srcAtlasVert(:,1)<0);
    AtlRInd=find(srcAtlasVert(:,1)>=0);
    % If hemispheres are connected: process all at once
    if any(isConnected)
        
        rHsrc  = 1:nSrc;
        rHdest = 1:nDest;
        lHsrc  = [];
        lHdest = [];
    end
    % Re-interpolate using the Squares and the shepards algorithm. First
    % find interpolation between subject and BrainSuiteAtlas1 and then find
    % interpolation matrix between BrainSuiteAtlas1 to Default Anatomy in
    % BrainStorm
    WmatTmpSrc2Atlas = bst_shepards(srcAtlasVert(AtlRInd,:), srcVert(rHsrc,:), nbNeighbors, 0);
    WmatTmpAtlas2Dest = bst_shepards(destVert(rHdest,:), destAtlasVert(AtlRInd,:), nbNeighbors, 0);
    WmatTmp=WmatTmpAtlas2Dest*WmatTmpSrc2Atlas;
    Wmat(rHdest,rHsrc) = WmatTmp;
    if ~isempty(lHdest)
        WmatTmpSrc2Atlas = bst_shepards(srcAtlasVert(AtlLInd,:), srcVert(lHsrc,:), nbNeighbors, 0);
        WmatTmpAtlas2Dest = bst_shepards(destVert(lHdest,:), destAtlasVert(AtlLInd,:), nbNeighbors, 0);
        WmatTmp=WmatTmpAtlas2Dest*WmatTmpSrc2Atlas;
        Wmat(lHdest,lHsrc) = WmatTmp;
    end
end


% ===== DEFAULT METHOD: BAD =====
% Else: Compute interpolation matrix
if isempty(Wmat)
    % Close all figures
    bst_memory('UnloadAll', 'Forced');
    % Warning: bad technique
    java_dialog('warning', ['This projection method you are about is outdated and inaccurate.' 10 10 ...
        'For accurate results, please consider using FreeSurfer for the MRI' 10 ...
        'segmentation, because it also generates registered spheres that we' 10 ...
        'can in Brainstorm for the the inter-subject co-registration.' 10 10 ...
        'More information in the FreeSurfer tutorial on the Brainstorm website: ' 10 ...
        'http://neuroimage.usc.edu/brainstorm/Tutorials/LabelFreeSurfer']);
    
    % === GET FIDUCIALS ===
    % Fiducials 3D positions are saved in the subject's MRI structure
    if ~isempty(sSrcSubj.Anatomy) && ~isempty(sDestSubj.Anatomy)
        % Get MRI filenames
        srcMriFile  = file_fullpath(sSrcSubj.Anatomy(sSrcSubj.iAnatomy).FileName);
        destMriFile = file_fullpath(sDestSubj.Anatomy(sDestSubj.iAnatomy).FileName);
        % Load NCS structures from MRIs (contains the fiducials AC,PC,IH)
        srcMri  = load(srcMriFile,  'NCS', 'SCS');
        destMri = load(destMriFile, 'NCS', 'SCS');
        % Check NCS and SCS fields
        isMissingSrc = ~isfield(srcMri,  'NCS') || ~all(isfield(srcMri.NCS,  {'AC','PC','IH'}))    || isempty(srcMri.NCS.AC)   || isempty(srcMri.NCS.PC)  || isempty(srcMri.NCS.IH) || ...
            ~isfield(srcMri,  'SCS') || ~all(isfield(srcMri.SCS,  {'NAS','LPA','RPA'})) || isempty(srcMri.SCS.NAS)  || isempty(srcMri.SCS.LPA) || isempty(srcMri.SCS.RPA);
        isMissingDest = ~isfield(destMri,  'NCS') || ~all(isfield(destMri.NCS,  {'AC','PC','IH'}))    || isempty(destMri.NCS.AC)   || isempty(destMri.NCS.PC)  || isempty(destMri.NCS.IH) || ...
            ~isfield(destMri,  'SCS') || ~all(isfield(destMri.SCS,  {'NAS','LPA','RPA'})) || isempty(destMri.SCS.NAS)  || isempty(destMri.SCS.LPA) || isempty(destMri.SCS.RPA);
        % If missing SOURCE fiducials
        if isMissingSrc && ~file_compare(srcMriFile, destMriFile)
            if isInteractive
                isStop = java_dialog('confirm', ['Warning: some fiducial points have not been marked on the source MRI.' 10 ...
                    'Without those fiducials, the projection will be much less accurate.' 10 10 ...
                    'Start MRI Viewer and define those points now ?'], 'Project sources');
                % Stop and start MRI Viewer
                if isStop
                    view_mri(srcMriFile, 'EditFiducials');
                    return
                end
            else
                bst_report('Error', 'process_project_sources', [], 'Some fiducial points have not been marked on the source MRI.');
                return
            end
        end
        % If missing DESTINATION fiducials
        if isMissingDest && ~file_compare(srcMriFile, destMriFile)
            if isInteractive
                isStop = java_dialog('confirm', ['Warning: some fiducial points have not been marked on the destination MRI.' 10 ...
                    'Without those fiducials, the projection will be much less accurate.' 10 10 ...
                    'Start MRI Viewer and define those points now ?'], 'Project sources');
                % Stop and start MRI Viewer
                if isStop
                    view_mri(destMriFile, 'EditFiducials');
                    return
                end
            else
                bst_report('Error', 'process_project_sources', [], 'Some fiducial points have not been marked on the destination MRI.');
                return
            end
        end
        % If everything is set correctly
        if ~isMissingSrc && ~isMissingDest
            srcSurfMat.NCS  = srcMri.NCS;
            srcSurfMat.SCS  = srcMri.SCS;
            destSurfMat.NCS = destMri.NCS;
            destSurfMat.SCS = destMri.SCS;
        end
    else
        errMsg = 'No available MRI for source or destination subject.';
        if isInteractive
            bst_error(errMsg, 'Project sources', 0);
        else
            bst_report('Error', 'process_project_sources', [], errMsg);
        end
        return
    end
    
    % === COMPUTE INTERPOLATION ===
    if isInteractive
        bst_progress('start', 'Project sources', 'Computing transformation...');
    end
    % Set the number of nearest neighbors to find
    %nbNeigh = 8 * ceil(length(destSurfMat.Vertices) / length(srcSurfMat.Vertices));
    nbNeigh = 8;
    % Interpolation cortex => cortex
    Wmat = tess_interp_cortex(destSurfMat, srcSurfMat, isSameSubject, nbNeigh);
end

% ===== SAVE INTERPOLATION =====
% Save interpolation in surface file, for future use
if isSaveInterp
    s.tess2tess_interp.Wmat      = Wmat;
    s.tess2tess_interp.Signature = Signature;
    bst_save(file_fullpath(srcSurfFile), s, 'v7', 1);
end


