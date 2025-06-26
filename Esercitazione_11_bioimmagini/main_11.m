%% caricamento dati

load 'MRI_es.mat';
imm = cell(1,2);
imm{1} = MRI_z90;
imm{2} = MRI_x100;
titoli.immagini = {'MRI_z90', 'MRI_x100'};

%% plot immagini e istogrammi

figure();
for i = 1:2
    subplot(2,2,i)
    imagesc(imm{i});
    title(titoli.immagini{i});
    colormap(gray)

    subplot(2,2,i+2)
    imhist(imm{i})
end

%% filtraggio immagini

titoli.filtri = {'filtro medio dimensione 3', 'filtro mediano dimensione 5', ...
    'filtro gradiente', 'filtro gaussiano'};

% applicazione filtri
imm_filt = cell(4, 2); % tipo filtro (medio 3, mediano 5, gauss, grad) x immagini (z90, x100)

for i = 1:2
    imm_filt{1, i} = imfilter(imm{i}, fspecial('average', 3)); % filtro medio
    imm_filt{2, i} = medfilt2(imm{i}, [5 5]); % filtro mediano
    imm_filt{3, i} = imfilter(imm{i}, fspecial('prewitt')); % grad y
    imm_filt{4, i} = imgaussfilt(imm{i}, 1, 'FilterSize', 3);
end 

%plot immagini filtrate

for i = 1:2
    figure
    sgtitle([titoli.immagini{i}, ' ', 'filtrata'])
    for j=1:4
        subplot(2,2,j);
        imagesc(imm_filt{j,i})
        colormap(gray)
        title(titoli.filtri{j})
    end
end


%% calcolo TDF

TDF = cell(3,2); % righe: originale, mediano 5, gradente. colonne: immagini

% calcolo della TDF 2D

for i = 1:2
    TDF{1,i} = fft2(imm{i});
    TDF{2,i} = fft2(imm_filt{2, i});
    TDF{3,i} = fft2(imm_filt{4, i});
end

%% plot TDF

% plot di moduli e fase

titoli.TDF = {'FDT immagini originali', 'FDT immagini filtrate (filtro mediano dim 5)', ...
    'FDT immagini filtrate (filtro gradiente)'};
N = zeros(4, 1);
[N(1) , N(2)] = size(imm{1});
[N(3) , N(4)] = size(imm{2});
freq= cell(2, 2); % righe = freq_y, freq_x. colonne = immagine

for i = 1:4
    freq{i} = linspace(-0.5, 0.5, N(i));
end

% plot 2D
for i = 1:3
    figure()
    sgtitle(titoli.TDF{i})

    for j =1:2
        modulo = abs(TDF{i,j});
        fase = angle(TDF{i,j});
        
        subplot(2,2,j)
        imagesc(fftshift(log(1 + modulo)));
        xticks(linspace(1, N(2*j), 5))
        xticklabels([-0.5, - 0.25,  0, 0.25, 0.5])
        yticks(linspace(1, N(2*j - 1), 5))
        yticklabels([-0.5, - 0.25,  0, 0.25, 0.5])
        xlabel("fx")
        ylabel("fy")
        title(['modulo', ' ', titoli.immagini{j}])
        colormap(gray)

        subplot(2,2,j+2)
        imagesc(fftshift(fase));
        xticks(linspace(1, N(2*j), 5))
        xticklabels([-0.5, - 0.25,  0, 0.25, 0.5])
        yticks(linspace(1, N(2*j - 1), 5))
        yticklabels([-0.5, - 0.25,  0, 0.25, 0.5])
        xlabel("fx")
        ylabel("fy")
        title(['fase', ' ', titoli.immagini{j}])
    end
end

% plot 3D
for i = 1:3
    figure()
    sgtitle(titoli.TDF{i})

    for j =1:2
        modulo = abs(TDF{i,j});
        fase = angle(TDF{i,j});
        
        subplot(2,2,j)
        surf(freq{2,j}, freq{1,j}, fftshift(log(1 + modulo)));
        title(['modulo', ' ', titoli.immagini{j}])
        xlabel("fx")
        ylabel("fy")
        shading interp
        
        subplot(2,2,j+2)
        surf(freq{2,j}, freq{1,j}, fftshift(fase));
        title(['fase', ' ', titoli.immagini{j}])
        xlabel("fx")
        ylabel("fy")
        shading interp
    end
end

%% segmentazione 

% ROI create con app image segmenter
ROI = cell(3, 2); %Sfondo, tessuto1, tessuto2
ROI{1,1} = imm{1}(1:30, 186:233);
ROI{2,1} = imm{1}(126:138, 143:162);
ROI{3,1} = imm{1}(37:48, 103:109);
ROI{1,2} = imm{2}(1:30, 186:233);
ROI{2,2} = imm{2}(140:159, 100:113);
ROI{3,2} = imm{2}(103:117, 159:165);

titoli.ROI = {'Sfondo', 'Tessuto1', 'Tessuto2'};

% plot ROI
figure()
for i = 1:2
    subplot(3,2,i)
    imagesc(ROI{1,i});
    title(['ROI ', titoli.ROI{1}, ' ', titoli.immagini{i}]);
    colormap(gray)

    subplot(3,2,i+2)
    imagesc(ROI{2,i});
    title(['ROI ', titoli.ROI{2}, ' ', titoli.immagini{i}]);
    colormap(gray)

    subplot(3,2,i+4)
    imagesc(ROI{3,i});
    title(['ROI ', titoli.ROI{3}, ' ', titoli.immagini{i}]);
    colormap(gray)
end

%% Thresholding
soglie = [50 75 100 130 160 175 190 220 250];
figure;

for i = 1:length(soglie)
    subplot(3, 3, i);
    img_thresh_emp = imm{1} > soglie(i);
    imshow(img_thresh_emp), 
    title(['Soglia = ', num2str(soglie(i)), ' ', titoli.immagini{1}]);
end

figure;
for i = 1:length(soglie)
    subplot(3, 3, i);
    img_thresh_emp = imm{2} > soglie(i);
    imshow(img_thresh_emp), 
    title(['Soglia = ', num2str(soglie(i)), ' ', titoli.immagini{2}]);
end

soglie = cell(1, 2);
soglie{1} = 175; % soglia per immagine 1
soglie{2} = 130; % soglia per immagine 2

%% SNR e CNR

SNR_originale = zeros(1, 2); % righe: immagine, colonne: ROI
CNR_originale = zeros(1, 2); % righe: immagine, colonne: ROI

for i = 1:2
    SNR_originale(i, 1) = mean(ROI{2, i}(:)) / std(double(ROI{1, i}(:))); % SNR Tessuto1
    CNR_originale(i, 1) = (mean(ROI{2, i}(:)) - mean(ROI{3, i}(:))) / std(double(ROI{1, i}(:))); % CNR Tessuto1
end

imm_thresh = cell(1, 2); % immagini sogliate
imm_thresh{1} = imm{1} > soglie{1};
imm_thresh{2} = imm{2} > soglie{2};

ROI_thresh = cell(3, 2); %Sfondo, tessuto1, tessuto2
ROI_thresh{1,1} = imm_thresh{1}(1:30, 186:233);
ROI_thresh{2,1} = imm_thresh{1}(126:138, 143:162);
ROI_thresh{3,1} = imm_thresh{1}(37:48, 103:109);
ROI_thresh{1,2} = imm_thresh{2}(1:30, 186:233);
ROI_thresh{2,2} = imm_thresh{2}(140:159, 100:113);
ROI_thresh{3,2} = imm_thresh{2}(103:117, 159:165);

SNR_thresh = zeros(1, 2); % righe: immagine, colonne: ROI
CNR_thresh = zeros(1, 2); % righe: immagine, colonne
for i = 1:2
    SNR_thresh(i, 1) = mean(ROI_thresh{2, i}(:)) / std(double(ROI_thresh{1, i}(:))); % SNR Tessuto1
    CNR_thresh(i, 1) = (mean(ROI_thresh{2, i}(:)) - mean(ROI_thresh{3, i}(:))) / std(double(ROI_thresh{1, i}(:))); % CNR Tessuto1
end

disp(titoli.immagini{1});
disp(['SNR originale: ', num2str(SNR_originale(1))]);
disp(['CNR originale: ', num2str(CNR_originale(1))]);
disp(['SNR sogliato: ', num2str(SNR_thresh(1))]);
disp(['CNR sogliato: ', num2str(CNR_thresh(1))]);
disp('');
disp(titoli.immagini{2});
disp(['SNR originale: ', num2str(SNR_originale(2))]);
disp(['CNR originale: ', num2str(CNR_originale(2))]);
disp(['SNR sogliato: ', num2str(SNR_thresh(2))]);
disp(['CNR sogliato: ', num2str(CNR_thresh(2))]);