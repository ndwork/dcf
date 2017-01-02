
function out = applyC2Grid_2D( F, kTraj, N, kCy, kCx, Cy, Cx )
  % out = applyC2Grid_2D( F, kTraj, N, kCx, kCy, Cx, Cy )
  %
  % Inputs:
  % F - 1D array of Fourier values
  % kTraj - Mx2 array where first/second column represents ky/kx location
  %   of trajectory
  % N - 2 element array specifying grid size [Ny Nx]
  %
  % Written by Nicholas Dwork - Copyright 2016
  % Based on codes written by John Pauly and Ethan Johnson

  kws = [ max(kCy)*N(1), max(kCx)*N(2) ];  % kernel widths
  rkws = round( kws );

  % Convert k-space coordinates to (non-integer) indexes into grid
  kTraj = mod( kTraj+0.5, 1 ) - 0.5;
  kTrajIndxs = zeros( size(kTraj) );
  kCoords = size2fftCoordinates( N );
  kTrajIndxs(:,1) = (1-ceil((N(1)+1)/2))/min(kCoords{1})*kTraj(:,1) + ceil((N(1)+1)/2);
  kTrajIndxs(:,2) = (1-ceil((N(2)+1)/2))/min(kCoords{2})*kTraj(:,2) + ceil((N(2)+1)/2);
  kCyIndxs = kCy * N(1);
  kCxIndxs = kCx * N(2);

  out = zeros(N);
  for dkx = -rkws(2):rkws(2)
    nearOutKxIndxs = round( kTrajIndxs(:,2) + dkx );
    kxDists = abs( kTrajIndxs(:,2) - nearOutKxIndxs );
    cValsKx = interp1( kCxIndxs, Cx, kxDists, 'linear', 0 );
    nearOutKxIndxs = mod( nearOutKxIndxs-1, N(2) ) + 1;

    for dky = -rkws(1):rkws(1)
      nearOutKyIndxs = round( kTrajIndxs(:,1) + dky );
      kyDists = abs( kTrajIndxs(:,1) - nearOutKyIndxs );
      cValsKy = interp1( kCyIndxs, Cy, kyDists, 'linear', 0 );
      nearOutKyIndxs = mod( nearOutKyIndxs-1, N(1) ) + 1;

      out = out + sparse( nearOutKyIndxs, nearOutKxIndxs, ...
       F .* cValsKy .* cValsKx, N(1), N(2) );
    end
  end

end
