
function [kTraj,iGridFVals,N,psfMask] = loadDataCase( datacase, varargin )
  % [kTraj,iGridFVals] = loadDataCase( datacase [, ds ] )
  %
  % ds is the amount to downsample the data

  dataDir = './data/';

  p = inputParser;
  p.addOptional( 'ds', 1, @isnumeric );
  p.parse(varargin{:});
  ds = p.Results.ds;

  psfMask = [];

  switch datacase
    case 1
      dataFile = 'data_2dspiralNav.mat';
      load( [dataDir dataFile] );
      sData = size( d );
      nCoils = sData(3);
      iGridFVals = zeros( prod(sData(1:2)), nCoils ); 
      for cIndx=1:nCoils
        iGridFVals(:,cIndx) = reshape( d(:,:,cIndx), [prod(sData(1:2)) 1] );
      end
      ks = k(:);
      kTraj = [ real(ks) imag(ks) ];

      N = [ 140 70 ];
  end

  %figure; scatter( kTraj(:,1), kTraj(:,2), 8, 'k', 'filled' );
end