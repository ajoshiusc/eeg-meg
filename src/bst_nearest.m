function [I, dist, dt] = bst_nearest(refVert, testVert, K, isProgress, dt)
% BST_NEAREST: K-Nearest neighbor search
% In this implementation, for each of the test vertices, we first find closest point 
% in the reference vertices. Then, a patch is expanded in the reference
% vertices till there are minimum 3K times vertices in the patch. On an average, 
% this patch contains 20K times vertices. The K nearest neighbors are found in this patch.
%
% USAGE:  I = bst_nearest(refVert, testVert, K=1, isProgress=1, dt=[])
%
% INPUT:
%    - refVert    : [Nx3] list of reference 3D points
%    - testVert   : [Mx3] list of 3D points for which we want the nearest vertex in refVert
%    - K          : Number of nearest neighbors we want, default=1
%    - isProgress : If 1, show a progress bar
%    - dt         : Delaunay triangulation returned by previous call
% OUTPUT:
%    - I : [MxK] search matrix
%    - d : [MxK] distance matrix

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
% Authors: Anand A. Joshi 2015

% Parse inputs
if (nargin < 5) || isempty(dt)
    dt = [];
end
if (nargin < 4) || isempty(isProgress)
    isProgress = 1;
end
if (nargin < 3) || isempty(K)
    K = 1;
end

if isProgress
    isPreviousBar = bst_progress('isVisible');
    if isPreviousBar
        initPos = bst_progress('get');
        bst_progress('text', 'Nearest neighbor search...');
    else
        bst_progress('start', 'Nearest neighbor', 'Nearest neighbor search...');
    end
end

if isProgress
    bst_progress('start', 'Nearest neighbor', 'Nearest neighbor search...', 1, 100);
end
 

Pts=[refVert;testVert];Pts=Pts+1000*eps*rand(size(Pts));
Ptsref=[refVert];Ptsref=Ptsref+1000*eps*rand(size(Ptsref));

fullt = delaunayTriangulation(Pts);
reft = delaunayTriangulation(Ptsref);

ed=fullt.edges;
Adj=sparse(ed(:,1),ed(:,2),1,length(fullt.Points),length(fullt.Points));
Adj=Adj+Adj';

ed=reft.edges;
Adjref=sparse(ed(:,1),ed(:,2),1,length(reft.Points),length(reft.Points));
Adjref=Adjref+Adjref';

Adjref2=0*Adj;
Adjref2(1:length(refVert),1:length(refVert))=Adjref+speye(size(Adjref));

AdjN=Adj(1:length(refVert),length(refVert)+1:end);
nbrs=full(sum(AdjN,1));Adj1=Adj+speye(size(Adj));
while min(nbrs)<3
    Adj=Adj1*Adj;Adj=(Adj>0);
    AdjN=Adj(1:length(refVert),length(refVert)+1:end);
    nbrs=full(sum(AdjN,1));
end
AdjN=Adj(1:length(refVert),length(refVert)+1:end);
nbrs=sum(AdjN,1);
iter=0;
while min(nbrs)<3*K
    Adj=Adjref2*Adj;Adj=(Adj>0);
    AdjN=Adj(1:length(refVert),length(refVert)+1:end);
    nbrs=sum(AdjN,1);%min(nbrs);
    iter=iter+1;
    if isProgress 
        pos = iter;
        bst_progress('set', pos);
    end
    
end
Adj=AdjN;
[r,c]=find(Adj);
Xref=0*Adj;
Ytst=0*Adj;

Xref=sparse(r,c,refVert(r,1),size(Xref,1),size(Xref,2));
Yref=sparse(r,c,refVert(r,2),size(Xref,1),size(Xref,2));

Xtst=sparse(r,c,testVert(c,1),size(Xref,1),size(Xref,2));
Ytst=sparse(r,c,testVert(c,2),size(Xref,1),size(Xref,2));

D2=(Xref-Xtst).^2 + (Yref-Ytst).^2;

D2inv=spfun(@(x)(1./x),D2);
%[B,I]=sort(D2inv,2);
I=zeros(size(testVert,1),K);dist=I;
nTest=size(D2inv,2);
pos=20;
for jj=1:nTest
    
    if isProgress && (jj/(nTest)*100 > pos)
        pos = ceil(jj/nTest*100);
        bst_progress('set', pos);
    end
    
    [dd,ii]=sort(D2inv(:,jj));
    I(jj,:)=ii(end-K+1:end);
    dist(jj,:)=sqrt(1./dd(end-K+1:end));
end

% Close progress bar
if isProgress && ~isPreviousBar
    bst_progress('stop');
end


