
function out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx, gridKs )
  % out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, gridKs ] )
  %
  % Written by Nicholas Dwork - Copyright 2016

  if nargin < 8
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};
    [gridKx,gridKy] = meshgrid(gridKx,gridKy);
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
  end

  nTraj = size(kTraj,1);
  kws = [ max(kCy), max(kCx) ];
  kDistThreshY = kws(1);
  kDistThreshX = kws(2);
  tmp = cell(nTraj,1);
  for trajIndx=1:nTraj
    distsKy = abs( kTraj(trajIndx,1) - gridKy );
    distsKx = abs( kTraj(trajIndx,2) - gridKx );
    shortDistIndxs = find( distsKy < kDistThreshY & ...
                           distsKx < kDistThreshX );
    shortDistsKy = distsKy( shortDistIndxs );
    shortDistsKx = distsKx( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    %out(shortDistIndxs) = out(shortDistIndxs) + ...
    %  F(trajIndx) * ( CValsY .* CValsX );

    tmp{trajIndx} = struct( ...
      'indxs', shortDistIndxs, ...
      'values', F(trajIndx) * ( CValsY .* CValsX ) ...
    );
  end

  out = zeros( size(gridKy) );
  for trajIndx=1:nTraj
    out( tmp{trajIndx}.indxs ) = out( tmp{trajIndx}.indxs ) + ...
      tmp{trajIndx}.values;
  end

  onesCol = ones(nTraj,1);
  for dim=1:2
    alt = zeros( size(kTraj) );
    alt(:,dim) = onesCol;

    for altDir=[-1 1]
      newTraj = kTraj + altDir*alt;
      if altDir < 0
        newTrajIndxs = find( newTraj(:,dim) > -0.5-kws(dim) );
      else
        newTrajIndxs = find( newTraj(:,dim) < 0.5+kws(dim) );
      end

      newTraj = newTraj( newTrajIndxs, : );
      for i=1:numel(newTrajIndxs)
        trajIndx = newTrajIndxs(i);
        NewDistsKy = abs( newTraj(i,1) - gridKy );
        NewDistsKx = abs( newTraj(i,2) - gridKx );
        NewShortDistIndxs = find( NewDistsKy < kDistThreshY & ...
                                  NewDistsKx < kDistThreshX );
        NewShortDistsKy = NewDistsKy( NewShortDistIndxs );
        NewShortDistsKx = NewDistsKx( NewShortDistIndxs );
        NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
        NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
        out(NewShortDistIndxs) = out(NewShortDistIndxs) + ...
          F(trajIndx) * ( NewCValsY .* NewCValsX );
      end
    end
  end

%   corners = { [-1 -1], [-1 1], [1 -1], [1 1] };
%   for cIndx = 1:numel(corners)
%     newTraj = kTraj + ones(size(kTraj)) .* repmat(corners, [nTraj 1]);
%   end
end
