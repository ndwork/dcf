
function makeResults
  close all;  drawnow;  clear;  rng(1);

  outDir = './results';
  algorithms = { 'FP', 'tbLSDC', 'rtbLSDC', 'gbLSDC' };
  goldAlg = 'FP';
  noiseScales = 0:0.01:0.1;

  datacases = makeDatacases();


  goldIndx = find( strcmp( algorithms, goldAlg ) );
  if goldIndx ~= 1
    tmp = algorithms{1};
    algorithms{1} = algorithms{goldIndx};
    algorithms{goldIndx} = tmp;
  end

  if ~exist(outDir,'dir'), mkdir( outDir ); end

  outFile = [outDir, '/makeResultsOut.csv'];
  outFileID = fopen( outFile, 'w' );
  fprintf( outFileID, 'Datacase, Algorithm, MSE, min weight, max weight\n' );

  trajDir = [outDir, '/trajectories'];
  if ~exist(trajDir,'dir'), mkdir( trajDir ); end;

  for caseIndx=1:numel(datacases)
    disp([ 'Working on case ', num2str( caseIndx ) ]);
    caseDir = [ outDir, '/case', num2str(caseIndx,'%3.3i'), '/' ];

    if isfield( datacases{caseIndx}, 'datacaseIndx' )

      [kTraj,iGridFVals,N,psfMask] = loadDataCase( ...
        datacases{caseIndx}.datacaseIndx, ...
        datacases{caseIndx}.downsample );

    else

      switch datacases{caseIndx}.trajType
        case 'poissonDisc'
          kTraj = makeTrajPts( 2, datacases{caseIndx}.trajType, ...
            datacases{caseIndx}.radius );
        case 'propeller'
          kTraj = makeTrajPts( 2, datacases{caseIndx}.trajType, ...
            datacases{caseIndx}.nReadout, datacases{caseIndx}.nLines, ...
            datacases{caseIndx}.dkLine, datacases{caseIndx}.nAngles );
        otherwise
          error('unrecognized trajectory type');
      end

      N = [128 128];
      radImg = makeRadialImg(2*N);
      %psfMask = radImg < abs(max(N));
      psfMask = ones(2*N);
    end

    figH = figure();
    scatter( kTraj(:,2), kTraj(:,1), 8, 'k', 'filled' );
    axis([-0.55 0.55 -0.55 0.55]);  axis('equal');
    set( gca, 'xtick', [], 'ytick', [] );
    img = frame2im(getframe(figH));
    trajImgFile = [trajDir,'/traj_case', num2str(caseIndx,'%3.3i'), '.png'];
    imwrite( img, trajImgFile );
    close( figH ); drawnow

    minWeight = 99d99;
    maxWeight = -minWeight;
    minDiff = minWeight;
    maxDiff = -minDiff;
    MSEs = cell( numel(algorithms), 1 );
    noiseImgFile = [caseDir, '/noiseSensitivity.png'];
    for algIndx = 1:numel(algorithms)
      thisAlg = algorithms{algIndx};
      disp([ 'Working on case ', num2str( caseIndx ), ' ', thisAlg ]);

      resultDir = [ caseDir, '/', thisAlg, '/'];
      if ~exist(resultDir,'dir'), mkdir( resultDir ); end;

      weightsMatFile = [resultDir,'/weights_',thisAlg,'.mat'];
      if exist( weightsMatFile, 'file' )
        load( weightsMatFile );
      else
        [weights,~,~] = makePrecompWeights_2D( kTraj, N, 'alg', thisAlg, ...
          'psfMask', psfMask );
        save( weightsMatFile, 'weights' );
      end

      if algIndx==1
        goldWeights = weights;
      else
        goldWeightsDiff = weights - goldWeights;
        if min(goldWeightsDiff(:)) < minDiff
          minDiff = min(goldWeightsDiff(:));
        end
        if max(goldWeightsDiff(:)) > maxDiff
          maxDiff = max(goldWeightsDiff(:));
        end
      end;

      if min(weights(:)) < minWeight
        minWeight = min( weights(:) );
      end
      if max(weights(:)) > maxWeight
        maxWeight = max( weights(:) );
      end

      if isempty( psfMask ), psfMask = ones(2*N); end;
      psf = grid_2D( ones(size(weights)), kTraj, 2*N, weights ) .* psfMask;
      psf = psf ./ max( psf(:) );
      psf_dB = 20*log10( psf ) .* psfMask;
      psf_dB(~isfinite(psf_dB)) = max( psf_dB(:) );
      psf_dB = scaleImg( psf_dB, [-100 0] );
      imwrite( psf_dB, [ resultDir, '/psf_dB.jpg' ] );
      b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
      mse = sum( abs( psfMask(:).*psf(:) - psfMask(:).*b(:) ).^2 ) ./ ...
        sum(psfMask(:));

      fprintf( outFileID, ['%4.4i, ', thisAlg, ', %12.10f, %12.10f, %12.10f \n'], ...
        caseIndx, mse, min(weights(:)), max(weights(:)) );

      if ~exist( noiseImgFile, 'file' )
        MSEs{algIndx} = findNoiseSensitivity( kTraj, weights, ...
          N, psfMask, noiseScales );
      end

      reconFile = [ resultDir, '/recon.png' ];
      if exist( 'iGridFVals', 'var' ) && ~exist( reconFile, 'file' )
        nCoils = size( iGridFVals, 2 );
        gridded = zeros( [ N nCoils] );
        for coilIndx=1:nCoils
          d = iGridFVals(:,coilIndx);  d = d(:);
          gridded_2D = grid_2D( d, kTraj, N, weights );
          gridded(:,:,coilIndx) = gridded_2D;
        end
        recon = sqrt( sum( gridded.*conj(gridded), 3 ) );
        imwrite( scaleImg(recon), reconFile );
      end
    end

    if ~exist( noiseImgFile, 'file' )
      tblsdcMSEs = MSEs{ strcmp( algorithms, 'tbLSDC' ) };
      rtblsdcMSEs = MSEs{ strcmp( algorithms, 'rtbLSDC' ) };
      noiseH = figure;
      plotnice( noiseScales, tblsdcMSEs, 'k--', 'lineWidth', 2 );
      hold on;
      plotnice( noiseScales, rtblsdcMSEs, 'b', 'lineWidth', 2 );
      ax = gca;
      set( ax, 'ytick', [min(ax.YTick) max(ax.YTick)] );
      legend( 'tbLSDC', 'rtbLSDC', 'Location', 'northwest' );
      saveas( noiseH, noiseImgFile );
      close; drawnow;
    end

    for algIndx = 1:numel(algorithms)
      thisAlg = algorithms{algIndx};
      resultDir = [ caseDir, '/', thisAlg, '/'];
      weightsMatFile = [resultDir,'/weights_',thisAlg,'.mat'];
      load( weightsMatFile );

      figH = figure;
      scatter( kTraj(:,1), kTraj(:,2), 4, weights, 'filled' );
      colormap('jet');  colorbar();
      axis([-0.55 0.55 -0.55 0.55]);  axis('equal');
      set( gca, 'xtick', [], 'ytick', [] );
      set( gca, 'fontsize', 16, 'LineWidth', 1.5 );
      img = frame2im(getframe(figH));
      weightsImgFile = [resultDir,'/weightsScaled.png'];
      imwrite( img, weightsImgFile );
      close;

      figH = figure;
      scatter( kTraj(:,1), kTraj(:,2), 4, weights, 'filled' );
      colormap('jet');  colorbar();  caxis([minWeight maxWeight]);
      axis([-0.55 0.55 -0.55 0.55]);  axis('equal');
      set( gca, 'xtick', [], 'ytick', [] );
      set( gca, 'fontsize', 16, 'LineWidth', 1.5 );
      img = frame2im(getframe(figH));
      weightsImgFile = [resultDir,'/weights.png'];
      imwrite( img, weightsImgFile );
      close;

      psf = grid_2D( ones(size(weights)), kTraj, 2*N, weights );
      psf = psf ./ max( psf(:) );
      if isempty(psfMask), psfMask = ones(size(psf)); end;
      psf_dB = 20*log10( psf ) .* psfMask;

      if algIndx==1
        goldWeights = weights;
        goldPSF_dB = psf_dB;
      else
        goldWeightsDiff = abs(weights - goldWeights);
        figH = figure;
        scatter( kTraj(:,1), kTraj(:,2), 4, goldWeightsDiff, 'filled' );
        %colormap('jet');  colorbar();  caxis([minDiff maxDiff]);
        colormap('jet');  colorbar();
        caxis([-1d-4 1d-4]);
        axis([-0.55 0.55 -0.55 0.55]);  axis('equal');
        set( gca, 'xtick', [], 'ytick', [] );
        set( gca, 'fontsize', 16, 'LineWidth', 1.5 );
        img = frame2im(getframe(figH));
        weightsImgFile = [ resultDir, '/weightsDiff.png' ];
        imwrite( img, weightsImgFile );
        close;

        psfDiff_dB = abs(psf_dB - goldPSF_dB);
        psfDiff_dB = scaleImg( psfDiff_dB, [0 30] );
        imwrite( psfDiff_dB, [resultDir, '/psfDiff_dB.jpg'] );
      end

      if exist( 'iGridFVals', 'var' )
        nCoils = size( iGridFVals, 2 );
        gridded = zeros( [ N nCoils] );
        for coilIndx=1:nCoils
          d = iGridFVals(:,coilIndx);  d = d(:);
          gridded_2D = grid_2D( d, kTraj, N, weights );
          gridded(:,:,coilIndx) = gridded_2D;
        end
        recon = sqrt( sum( gridded.*conj(gridded), 3 ) );
        reconFile = [resultDir,'/recon.png'];
        if algIndx==1
          scaleMin = mean(recon(:)) - 2.5*std(recon(:));
          scaleMax = mean(recon(:)) + 2.5*std(recon(:));
          scaleMinMax = [scaleMin scaleMax];
        end
        imwrite( scaleImg(recon), reconFile );
        scaledRecon = scaleImg( recon, scaleMinMax );
        scaledReconFile = [resultDir,'/scaledRecon.png'];
        imwrite( scaledRecon, scaledReconFile );
      end
    end

    clear iGridFVals
  end

  fclose( outFileID );
end


function MSEs = findNoiseSensitivity( kTraj, weights, N, psfMask, noiseScales )
  nSample = 10;

  nNoises = numel(noiseScales);
  MSEs = zeros( nNoises, 1 );
  nWeights = numel(weights);

  idealPSF = zeros(2*N);  idealPSF(1)=1;
  idealPSF=fftshift(idealPSF);
  idealPSF = idealPSF .* psfMask;
  for noiseIndx = 1:nNoises
    disp([ 'Working on noise index ', num2str(noiseIndx), ' of ', ...
      num2str(nNoises) ]);
    theseMSEs = cell( nNoises, 1 );
    if isempty( psfMask ), psfMask = ones(size(psf)); end;

    parfor sampIndx=1:nSample
      noiseScale = noiseScales(noiseIndx);
      noise = randn(nWeights,1) * noiseScale;
      noisyFs = ones(size(weights)) + noise;

      psf = grid_2D( noisyFs, kTraj, 2*N, weights );
      psf = psf ./ max( psf(:) ) .* psfMask;
      theseMSEs{sampIndx} = norm( psf - idealPSF, 2 ) / nWeights;
    end
    theseMSEs = cell2mat( theseMSEs );

    MSEs(noiseIndx) = mean( theseMSEs );
  end

end



function datacases = makeDatacases()

  datacases{1} = struct( ...
    'datacaseIndx', 1, ...
    'downsample', 1 ...
  );

end

