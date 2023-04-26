clear;
% code for paper
% Arturo Canales-Benavides, Yue Zhuo, Andrea M. Amitrano, Minsoo Kim, Raul I. Hernandez-Aranda, P. Scott Carney, and Martin Schnell, "Accessible quantitative phase imaging in confocal microscopy with sinusoidal-phase synthetic optical holography," Appl. Opt. 58, A55-A64 (2019)

% Reconstruction parameters
ky = 311.3;  kx = 0.1; 
filter_x = 1/1; filter_y = 1/12;
termweight = 0.9; 

% Read hologram from file
mat_holo = imread('hologram.tif'); [Ny, Nx] = size(mat_holo);

% Calculate window function
[X,Y] = meshgrid(1:Nx, 1:Ny);
mat_window_x1 = 0.2 + 0.8 * cos( (X-Nx/2-1*(kx+1)) / (Nx*filter_x) * pi/2);
mat_window_x2 = 0.2 + 0.8 * cos( (X-Nx/2-2*(kx+1)) / (Nx*filter_x) * pi/2);
mat_window_y1 = 0.2 + 0.8 * cos( (Y-Ny/2-1*(ky+1)) / (Ny*filter_y) * pi/2);
mat_window_y2 = 0.2 + 0.8 * cos( (Y-Ny/2-2*(ky+1)) / (Ny*filter_y) * pi/2);
mat_mask_x1 = abs(X-Nx/2-1*(kx+1)) < Nx*filter_x;
mat_mask_x2 = abs(X-Nx/2-2*(kx+1)) < Nx*filter_x;
mat_mask_y1 = abs(Y-Ny/2-1*(ky+1)) < Ny*filter_y;
mat_mask_y2 = abs(Y-Ny/2-2*(ky+1)) < Ny*filter_y;
mat_window1 = mat_window_x1 .* mat_window_y1 .* mat_mask_x1 .* mat_mask_y1;
mat_window2 = mat_window_x2 .* mat_window_y2 .* mat_mask_x2 .* mat_mask_y2;

% Perform Fourier filtering
mat_fft = fftshift(fft2(mat_holo)); 
mat_fff1 = mat_fft .* mat_window1;
mat_fff2 = mat_fft .* mat_window2;
holo_reco1 = ifft2( fftshift(mat_fff1) );
holo_reco2 = ifft2( fftshift(mat_fff2) );

% Shift terms 1 and 2 to center
mat_phase1 = exp( -2i*pi*Y/Ny*1*ky) .* exp( -2i*pi*X/Nx*1*kx ) ;
mat_phase2 = exp( -2i*pi*Y/Ny*2*ky) .* exp( -2i*pi*X/Nx*2*kx ) ;
holo_reco1 = holo_reco1 .* mat_phase1;
holo_reco2 = holo_reco2 .* mat_phase2;

% Determination of phi_0
[histcnt,histbin] = histcounts( angle(holo_reco1(:)), 100);
[histmax,histmaxind] = max(histcnt);
phase_offset = histbin(histmaxind);
holo_reco1 = holo_reco1 .* exp(-1i * 1*phase_offset);
holo_reco2 = holo_reco2 .* exp(-1i * 2*phase_offset);

% Plot intermediate results
figure(2); imagesc(abs(mat_fft)); caxis( [ 0 max(max(abs(mat_fft)))/1000] ); colormap 'jet'; title 'FFT';
figure(3); hmax = max( max(max(abs(holo_reco1))), max(max(abs(holo_reco2))));
subplot(2,4,1); imagesc( abs(holo_reco1) ); caxis([ 0 hmax] ); title 'Term1 abs';
subplot(2,4,2); hist( abs(holo_reco1(:)), 100); title 'Hist Term1 abs';
subplot(2,4,3); imagesc( angle(holo_reco1) ); caxis([ -3.14 3.14]); title 'Term1 arg';
subplot(2,4,4); hist( angle(holo_reco1(:)), 100); title 'Hist Term1 arg';
subplot(2,4,5); imagesc( abs(holo_reco2) ); caxis([ 0 hmax] ); title 'Term2 abs';
subplot(2,4,6); hist( abs(holo_reco2(:)), 100); title 'Hist Term2 abs';
subplot(2,4,7); imagesc( angle(holo_reco2) ); caxis([ -3.14 3.14]); title 'Term2 arg';
subplot(2,4,8); hist( angle(holo_reco2(:)), 100); title 'Hist Term2 arg';
% Note: Adjust kx, ky until any gradient in the phase images of term 1 and
% 2 is eliminated. This can be recognized with the histogram of the phase
% that should only show counts at 0 and -pi/pi when the correct values for
% ky and kx were found (see screenshot).

% Threshold phase
mat_sign = ones(Ny,Nx); mat_sign( abs(angle(holo_reco1)) > pi/2 ) = -1; 
holo_reco1 = abs(holo_reco1) .* mat_sign;
mat_sign = ones(Ny,Nx); mat_sign( abs(angle(holo_reco2)) > pi/2 ) = -1; 
holo_reco2 = abs(holo_reco2) .* mat_sign;

% Final reconstruction
holo_reco = 1*real(holo_reco1) + termweight*1i* real(holo_reco2);
figure(1); imagesc(angle(holo_reco));

% Additionally, phase gradients owing to sample tilt can be corrected
phase_kx = -0.2; phase_ky = -0.2; phase_offs = -3; % manual adjustment needed
mat_phase_corr = exp( -2i*pi*Y/Ny*phase_ky) .* exp( -2i*pi*X/Nx*phase_kx ) * exp(1i * phase_offs)  ;
holo_reco = holo_reco .* mat_phase_corr;
holo_reco = conj(holo_reco); %optionally, invert phase contrast
figure(1); imagesc(angle(holo_reco));

