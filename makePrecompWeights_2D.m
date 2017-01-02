
function [weights,flag,res] = makePrecompWeights_2D( ...
  kTraj, N, varargin )
  % [weights,flag,res] = makePrecompWeights_2D( kTraj, N, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg, , ...
  %     'psfMask', psfMask ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   kTraj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column is the ky/kx location.
  %     The units are normalized to [-0.5,0.5)
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %
  % Optional Inputs:
  %   alpha - the oversampling factor > 1
  %   W - the window width in pixels
  %   nC - the number of points to sample the convolution kernel
  %   alg - a string specifying the algorithm to use
  %     FP - specifies fixed point iteration
  %     gbLSDC - minimizes in the frequency domain
  %     rtbLSDC (default) - robust least squares
  %     tbLSDC - specifies least squares
  %     WLS - minimizes in the space domain
  %   nIter - specifies the number of iterations of fp method
  %   psfMask - only used by space domain optimizations
  %
  % Outputs:
  %   weights - 1D array with density compensation weights
  %
  % Optional Outputs:
  %   flag - flag describing results of optimization (see lsqr
  %     documentation)
  %   res - residual
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultAlg = 'rtbLSDC';
  defaultNIter = 15;
  radImg = makeRadialImg(2*N);
	defaultPsfMask = radImg < abs(min(N));
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.addParameter( 'alg', defaultAlg );
  p.addParameter( 'nIter', defaultNIter, checknum );
  p.addParameter( 'psfMask', defaultPsfMask );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  alg = p.Results.alg;
  nIter = p.Results.nIter;
  psfMask = p.Results.psfMask;

  flag = 0;
  res = 0;
  switch alg
    case 'gbLSDC'
      % Least squares in frequency domain on grid points
      [weights,flag,res] = makePrecompWeights_2D_gbLSDC( ...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'FP'
      % Pipe's method
      [weights,flag,res] = makePrecompWeights_2D_FP( ...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC, 'nIter', nIter );
      
    case 'tbLSDC'
      % Optimization analogy of Pipe's Method
      [weights,flag,res] = makePrecompWeights_2D_tbLSDC( ...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'rtbLSDC'
      % tbLSDC with non-negativity constraint
      [weights,flag,res] = makePrecompWeights_2D_rtbLSDC( ...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'WLS'
      [weights,flag,res] = makePrecompWeights_2D_WLS( ...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC, 'psfMask', psfMask );
      
    otherwise
      error('makePrecompWeights: Algorithm not recognized');
  end

end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_gbLSDC( ...
  traj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;

  % Make the Kaiser Bessel convolution kernel
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  %cGrid = 2*N;
  cGrid = nGrid;

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, cGrid );
      out = applyCT_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
    else
      out = applyC_2D( in, traj, cGrid, kCy, kCx, Cy, Cx );
      out = out(:);

      disp(['makePrecompWeights_2D_gbLSDC working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b = ones(cGrid);
  tolerance = 1d-6;
  maxIter = 1000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), ...
    tolerance, maxIter );

  radialImg = makeRadialImg( 2*N );
  mask = double( radialImg <= min(N) );
  scale = showPSF( weights, traj, N, mask, 'imgTitle', 'gbLSDC');
  weights = scale * weights;
  close;
end


function [weights,flag,res] = makePrecompWeights_2D_FP( ...
  traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  defaultNIter = 20;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.addParameter( 'nIter', defaultNIter, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  nIter = p.Results.nIter;

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  % Make the Kaiser Bessel convolution kernel
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  nTraj = size( traj, 1 );
  weights = ones( nTraj, 1 );

  flag = 1;
  for iteration=1:nIter
    if mod( iteration, 5 ) == 0,
      disp(['makePrecompWeights_2D_FP working on iteration ', ...
        num2str(iteration) ]);
    end

    oldWeights = weights;
    denom = applyC_2D( oldWeights, traj, N, kCy, kCx, Cy, Cx, traj );
    weights = oldWeights ./ denom;
  end

  radialImg = makeRadialImg( 2*N );
  mask = double( radialImg <= min(N) );

  scale = showPSF( weights, traj, N, mask, 'imgTitle', 'FP' );
  weights = scale * weights;
  close

  if nargout > 1
    flag = 0;
  end
  if nargout > 2
    res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );
  end
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_tbLSDC( ...
  traj, N, varargin )
  % Optimization analogy of Pipe's algorithm

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;

  % Make the Kaiser Bessel convolution kernel
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      out = applyCT_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_2D_tbLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  b = ones(size(traj,1),1);
  tolerance = 1d-5;
  maxIter = 5000;
  [weights,lsFlag,lsRes,lsIter,lsResVec,lsVec] = lsqr( @applyA, b(:), ...
    tolerance, maxIter );

  scale = showPSF( weights, traj, N, 'imgTitle', 'tbLSDC');
  weights = scale * weights;
  close;
end


function [weights,flag,residual] = makePrecompWeights_2D_rtbLSDC( ...
  traj, N, varargin )
  % Method of tbLSDC with a non-negativity constraint

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  iteration = 0;

  % Make the Kaiser Bessel convolution kernel
  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );
  Ny=N(1);  Nx=N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, trueAlpha, W, nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, trueAlpha, W, nC );

  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      out = applyCT_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if mod( iteration, 5 ) == 0,
        disp(['makePrecompWeights_2D_rtbLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  nTraj = size(traj,1);
  b = ones(nTraj,1);
  function y = linop_4_tfocs( in, mode )
    switch mode,
      case 0
        y = [numel(b), nTraj];
      case 1
        y = applyA( in );
      case 2
        y = applyA( in, 'transp' );
    end
  end
  %varargout = linop_test( @linop_4_tfocs );

  opts.alg = 'N83';
  opts = tfocs;
  opts.maxIts = 10000;
  opts.printEvery = 1;
  opts.tol = 1d-5;
  x0 = ones( nTraj, 1 );
  [weights,tfocsOptOut] = tfocs( smooth_quad, { @linop_4_tfocs, -b }, ...
    proj_Rplus, x0, opts );
  %weights = tfocs( smooth_huber(0.02), { @linop_4_tfocs, -b }, [], x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end
  if nargout > 2
    psf = applyA( weights );
    residual = norm( psf(:) - b(:), 2 ) / norm(b(:),2);
  end
  
tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  scale = showPSF( weights, traj, N, 'imgTitle', 'rtbLSDC');
  weights = scale * weights;
  close;
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_WLS( ...
  kTraj, N, varargin )

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  nGrid = ceil( 2 * N );
  radImg = makeRadialImg(nGrid);
	defaultPsfMask = radImg < abs(min(N));
  checknum = @(x) isnumeric(x) && isscalar(x) && (x > 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.addParameter( 'psfMask', defaultPsfMask );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  psfMask = p.Results.psfMask;

  iteration = 0;
  mask = psfMask;
  b=zeros(nGrid);  b(1,1)=1;  b=fftshift(b);
  mask(b>0) = sqrt( sum(mask(:)) );

  function out = applyWA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      in = reshape( in, nGrid );
      masked = mask .* in;
      out = iGrid_2D( masked, kTraj, 'alpha', alpha, 'W', W, 'nC', nC );
      out = real(out);
    else
      out = iGridT_2D( in, kTraj, nGrid, 'alpha', alpha, 'W', W, 'nC', nC );
      out = mask .* out;
      out = out(:);

      disp(['makePrecompWeights_2D_WLS working on iteration ', num2str(iteration) ]);
      iteration = iteration + 1;
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

tmp=0; tmp2=0; tmp3=0; tmp4=0; x0=0;

  tolerance = 1d-6;
  maxIter = 5000;
  [weights,lsFlag,lsRes] = lsqr( @applyWA, b(:), tolerance, maxIter );

  scale = showPSF( weights, kTraj, N, mask, 'imgTitle', 'WLS');
  weights = scale * weights;
  close;
end


