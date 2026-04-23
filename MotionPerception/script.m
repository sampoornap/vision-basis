% L1: Make of movie of two drifting gratings (sf: 2 c/deg, speed: 3 deg/sec, ori: –45 and +45 deg)
% L2: Make a movie of a plaid composed of the sum of both drifting gratings
% L3: Make a movie of a plaid whose components differ in orientation by an arbitrary amount
% L4: Make a movie of a plaid whose components produce transparent instead of coherent motion
% L5: Make a movie of 3-D luminance noise
% L6: Make a velocity filter in the 3–D Fourier domain and create filtered noise


% L1 

dims = [100 100];
pixSizeDeg = 0.05;
sf = 2; % cycles/deg
speed = 3; % deg/sec
ori1 = -45;
ori2 = 45;

nFrames = 60;
dt = 0.03;

movie1 = zeros([dims nFrames]);
movie2 = zeros([dims nFrames]);

for t = 1:nFrames
    phase = 2 * pi * sf * speed * dt * t;

    movie1(:,:,t) = giveSine(dims, pixSizeDeg, sf, ori1, phase);
    movie2(:,:,t) = giveSine(dims, pixSizeDeg, sf, ori2, phase);
end

figure;
for t = 1:nFrames
    imagesc(movie1(:,:,t));
    colormap gray;
    axis image off;
    title(['Frame ' num2str(t)]);
    drawnow;
end

%compare 2 gratings
figure;
for t = 1:nFrames
    subplot(1,2,1);
    imagesc(movie1(:,:,t), [0 1]);
    title('Grating 1');
    axis image off;

    subplot(1,2,2);
    imagesc(movie2(:,:,t), [0 1]);
    title('Grating 2');
    axis image off;

    colormap gray;
    drawnow;
end

% L2 

plaid = (movie1 + movie2) / 2;

% L3

ori1 = 0;
ori2 = 90; 

% L4

for t = 1:nFrames
    phase1 = 2 * pi * sf * speed * dt * t;
    phase2 = 2 * pi * sf * (speed * 0.5) * dt * t; % different speed

    g1 = giveSine(dims, pixSizeDeg, sf, ori1, phase1);
    g2 = giveSine(dims, pixSizeDeg, sf, ori2, phase2);

    transparent(:,:,t) = (g1 + g2) / 2;
end

% display plaid 
figure;
for t = 1:nFrames
    imagesc(plaid(:,:,t), [0 1]);
    colormap gray;
    axis image off;
    title('Plaid Motion');
    drawnow;
end

% L5

noise3D = randn([dims nFrames]);

%display noise
figure;
for t = 1:nFrames
    imagesc(noise3D(:,:,t));
    colormap gray;
    axis image off;
    title('3D Noise');
    drawnow;
end

%L6 

F = fftn(noise3D);

% first we create frequency grid
[fx, fy, ft] = ndgrid(linspace(-1,1,dims(1)), linspace(-1,1,dims(2)), linspace(-1,1,nFrames));

% velocity plane
vx = 0.2;
vy = 0.2;

% vel filter
H = exp(-((fx - vx*ft).^2 + (fy - vy*ft).^2) / 0.01);

F_filtered = F .* fftshift(H);

% back to space-time
filteredNoise = real(ifftn(F_filtered));

figure;
for t = 1:nFrames
    imagesc(filteredNoise(:,:,t));
    colormap gray;
    axis image off;
    title('Velocity-filtered noise');
    drawnow;
end