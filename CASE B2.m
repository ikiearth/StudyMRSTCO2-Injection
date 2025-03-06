mrstModule add ad-core compositional linearsolvers ad-props mrst-gui
close all;

% Mengaktifkan Gravitasi
gravity reset on

% === Fase 1: Geometri Reservoir & Petrofisika === %

% Dimensi Grid untuk Reservoir
I = 30; % Jumlah grid dalam arah I
J = 20; % Jumlah grid dalam arah J
K = 10; % Jumlah grid dalam arah K

gridDimension = [I, J, K];

pdim = [2000, 1000, 199];

% Membuat Grid Cartesian
G = cartGrid(gridDimension, pdim); % Membuat grid Cartesian
G = computeGeometry(G);            % Menghitung geometri grid

% Permeabilitas Reservoir
perm_mean_value = 13;  % Permeabilitas rata-rata (mD)
perm_std_value = 12;   % Standar deviasi permeabilitas (mD)
perm_min_value = 11;   % Permeabilitas minimum (mD)
perm_max_value = 15;   % Permeabilitas maksimum (mD)

% Konversi ke satuan MRST (meter persegi)
perm_mean = perm_mean_value * milli * darcy; % Permeabilitas rata-rata dalam m^2
perm_std = perm_std_value * milli * darcy;   % Standar deviasi dalam m^2

% Membuat distribusi normal untuk permeabilitas
perm_horizontal = perm_mean + perm_std * randn(G.cells.num, 1);

% Membatasi nilai permeabilitas horizontal
perm_horizontal = max(min(perm_horizontal, perm_max_value * milli * darcy), ...
                      perm_min_value * milli * darcy);

% Definisikan permeabilitas anisotropik
perm_vertical = 0.1 * perm_horizontal;                    % Permeabilitas vertikal 10% dari horizontal
perm = [perm_horizontal, perm_horizontal, perm_vertical]; % [Kx, Ky, Kz]

% Porositas Reservoir
poro_mean_value = 0.15; % Mean porositas
poro_std_value = 0.05;  % Standar deviasi porositas
poro_max_value = 0.16;  % Porositas maksimum
poro_min_value = 0.13;  % Porositas minimum

% Membuat porositas dengan distribusi normal
poro_mean = poro_mean_value;                         % Mean porositas
poro_std = poro_std_value;                           % Standar deviasi porositas
poro = poro_mean + poro_std * randn(G.cells.num, 1); % Porositas dengan distribusi normal

% Membatasi nilai porositas
poro = max(min(poro, poro_max_value), poro_min_value);

% Membuat Properti Batuan
rock = makeRock(G, perm, poro);
%% ------------------------------------------------------------------------------------------------------------------------

% === Phase 2: Initial State === %

% Inisialisasi Parameter Waktu Simulasi
simTimeProd = 3 * year;    % Total waktu simulasi produksi (5 tahun)
simTimeInj = 2 *year;      % Total waktu simulasi injeksi (5 tahun)
nstep   = 2 * day;         % Jumlah langkah waktu utama
refine  = 20;              % Refinement untuk langkah awal (langkah kecil di awal)

% Menghitung Volume Pori dan Laju Injeksi
pv      = poreVolume(G, rock);       % Menghitung volume pori dari grid dan batuan
injRate = 1 * sum(pv) / simTimeInj;  % Laju injeksi: fraksi volume pori dibagi total waktu simulasi

% Define flow fluid model
flowfluid = initSimpleADIFluid('phases', 'WOG', ...                % Phase Fluid
                           'mu',  [1, 10, 0.0001]*centi*poise, ... % Viskositas [air, minyak, gas]
                           'rho', [1000, 800, 10],...              % Densitas [air, minyak, gas]
                           'n', [2, 2, 2]);                        % Corey exponent

% Define connate saturations and residual saturations
Siw = 0.15;     % Connate water saturation
Sorw = 0.15;    % Residual oil saturation (water-oil system)
Sorg = 0.1;     % Residual oil saturation (gas-oil system)
Sgc = 0.05;     % Critical gas saturation

% Corey exponents and end-point relative permeabilities
krwro = 0.15;   % End-point kr of water at residual oil
krocw = 1.0;    % End-point kr of oil at connate water
krgro = 1.0;    % End-point kr of gas at residual oil
krogr = 1.0;    % End-point kr of oil at connate gas

% Define new relative permeability functions
flowfluid.krW = @(s) krwro * max((s - Siw) ./ (1 - Siw - Sorw), 0).^2;
flowfluid.krO_wo = @(s) krocw * max((1 - s - Sorw) ./ (1 - Siw - Sorw), 0).^2;

flowfluid.krO_go = @(s) krogr * max((1 - s - Sorg) ./ (1 - Siw - Sorg), 0).^2;
flowfluid.krG = @(s) krgro * max((s - Sgc) ./ (1 - Sgc - Siw - Sorg), 0).^2;

% Generate saturation range
Sw = linspace(Siw, 1 - Sorw, 1001);       % Water saturation
Sg = linspace(Sgc, 1 - Siw - Sorg, 1001); % Gas saturation
So = 1 - Sg - Sorg;                       % Oil saturation in gas-oil system

% Plot Relative Permeability Curve
figure('Position', [100, 100, 1000, 400]);

% Subplot kiri: Water-Oil relative permeability
subplot(1,2,1);
plot(Sw, flowfluid.krW(Sw), 'b', 'linewidth', 1.5); hold on;
plot(So, flowfluid.krO_wo(So), 'g', 'linewidth', 1.5);

title('Relative Permeability Curves: Water-Oil System');
xlabel('Water Saturation (S_w)');
ylabel('Relative Permeability (k_r)');
legend('Water', 'Oil');
grid on;
hold off;

% Subplot kanan: Gas-Oil relative permeability
subplot(1,2,2);
plot(Sg, flowfluid.krG(Sg), 'r', 'linewidth', 1.5); hold on;
plot(So, flowfluid.krO_go(So), 'g', 'linewidth', 1.5);

title('Relative Permeability Curves: Gas-Oil System');
xlabel('Gas Saturation (S_g)');
ylabel('Relative Permeability (k_r)');
legend('Gas', 'Oil');
grid on;
hold off;

% Compressibility
c_oil   = 5e-4 / barsa; % Compressibility minyak
c_water = 5e-6 / barsa; % Compressibility air
c_gas   = 3e-3 / barsa; % Compressibility gas
p_ref   = 275 * barsa;

% Formation volume factors
flowfluid.bO = @(p) exp((p - p_ref) * c_oil);   % Oil
flowfluid.bW = @(p) exp((p - p_ref) * c_water); % Water
flowfluid.bG = @(p) exp((p - p_ref) * c_gas);   % Gas

% Plot formation volume factors for oil (bO), water (bW), and gas (bG)
figure;

% Plot oil formation volume factor
p0 = linspace(100, 500, 10) * barsa;
plot(p0 / barsa, flowfluid.bO(p0), 'g-', 'DisplayName', 'b_O (Oil)')
hold on

% Plot water formation volume factor
plot(p0 / barsa, flowfluid.bW(p0), 'b--', 'DisplayName', 'b_W (Water)')

% Plot gas formation volume factor
plot(p0 / barsa, flowfluid.bG(p0), 'r-.', 'DisplayName', 'b_G (Gas)')

% Formatting the plot
xlabel('Pressure (bar)')
ylabel('Ratio')
title('Formation Volume Factors for Oil (bO), Water (bW), and Gas (bG)')
legend('show')
hold off

% Reference capillary entry pressure
pe_Wo = 0.5 * barsa; % Water-Oil entry capillary pressure
pe_Go = 0.3 * barsa; % Gas-Oil entry capillary pressure

% Capillary Pressure Functions Model Brooks-Corey
flowfluid.pcWO = @(sw) pe_Wo * ((max(min(sw - Siw, 1 - Siw - Sorw), 1e-5) / (1 - Siw - Sorw)).^-0.5);
flowfluid.pcGO = @(sg) pe_Go * ((max(min(sg - Sgc, 1 - Sgc - Siw - Sorg), 1e-5) / (1 - Sgc - Siw - Sorg)).^-0.5);

pcwo = flowfluid.pcWO(Sw);
pcgo = flowfluid.pcGO(Sg);

% Plot Capillary Pressure Curve
figure;

% Water-Oil
subplot(2,1,1);
plot(Sw, pcwo / barsa, 'b-', 'LineWidth', 2);
xlabel('Water Saturation (Sw)');
ylabel('Pc_{wo} [bar]');
title('Water-Oil Capillary Pressure');
axis([0 1 -5 30]);
grid on;

% Add annotations for residual saturations (Water-Oil)
yline(0, 'k--', 'LineWidth', 1); % Line at Pc = 0
xline(Siw, 'k:', 'LineWidth', 1); % Line for connate water saturation (Siw)
xline(1 - Sorw, 'r:', 'LineWidth', 1); % Line for residual oil saturation (Sorw)

% Add text annotations
text(Siw - 0.02, -3, '$S_{iw}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight','bold');
text(1 - Sorw - 0.02, -3, '$1 - S_{orw}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight','bold');

% Gas-Oil
subplot(2,1,2);
plot(Sg, pcgo / barsa, 'r-', 'LineWidth', 2);
xlabel('Gas Saturation (Sg)');
ylabel('Pc_{go} [bar]');
title('Gas-Oil Capillary Pressure');
axis([0 1 -5 30])
grid on;

% Add annotations for residual saturations (Gas-Oil)
yline(0, 'k--', 'LineWidth', 1); % Line at Pc = 0
xline(Sgc, 'k:', 'LineWidth', 1); % Line for critical gas saturation (Sgc)
xline(1 - Siw - Sorg, 'r:', 'LineWidth', 1); % Line for residual oil saturation (Sorg)

% Add text annotations
text(Sgc - 0.02, -3, '$S_{gc}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight','bold');
text(1 - Siw - Sorg - 0.05, -3, '$1 - S_{iw} - S_{org}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight','bold');
%% ------------------------------------------------------------------------------------------------------------------------

% === Phase 3: Set Up Fluid Model === %

% Menyiapkan Model Fluida
% Memuat model fluida komposisional
[fluid, info] = getBenchmarkMixture('lumped_1');

% Definisikan model EOS (Pengaturan EOS)
eosname = 'prcorr';
eosModel = EquationOfStateModel(G, fluid, eosname);

% Distribusi tekanan awal karena buoyancy
z = G.cells.centroids(:, 3);  % Koordinat vertikal
g = norm(gravity());          % Konstanta gravitasi

model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', true);

% Inisialisasi tekanan dan saturasi awal
% Set saturasi awal
s0 = [0.2, 0.5, 0.3]; % [sW, sO, sG]

% Definisikan kedalaman kontak
WOC_depth = 200; % meter, kedalaman Water-Oil Contact
GOC_depth = 50; % meter, kedalaman Gas-Oil Contact

% Definisikan tekanan di kedalaman referensi
ref_depth = 150;              % meter, kedalaman referensi untuk tekanan
ref_pressure = info.pressure; % tekanan di kedalaman referensi

% Inisialisasi saturasi
sW = zeros(size(z));
sO = zeros(size(z));
sG = zeros(size(z));

% Atur saturasi berdasarkan kedalaman
for i = 1:G.cells.num
    z = G.cells.centroids(i, 3);
    if z < GOC_depth
        % Di atas GOC: hanya gas
        state0.s(i, :) = [0, 0, 1];
    elseif z >= GOC_depth && z < WOC_depth
        % Antara GOC dan WOC: campuran sesuai s0
        state0.s(i, :) = s0;
    else
        % Di bawah WOC: hanya air
        state0.s(i, :) = [1, 0, 0];
    end
end

z = G.cells.centroids(:,3);

state0.pressure = ref_pressure + flowfluid.rhoOS * g * (z - ref_depth); % Tekanan minyak
state0 = initCompositionalState(G, state0.pressure , info.temp, s0, info.initial, model.EOSModel);

saturationWater = state0.s(:, 1); % Kolom 1 mewakili water
saturationOil = state0.s(:, 2);   % Kolom 1 mewakili minyak
saturationGas = state0.s(:, 3);   % Kolom 2 mewakili gas
saturationSum = saturationGas + saturationOil + saturationWater;

if any(abs(saturationSum - 1) > 1e-6)
    error('Saturasi tidak menjumlahkan hingga 1. Periksa inisialisasi!');
end

% Bulk Volume (V_b) dari grid
V_b = G.cells.volumes; % Volume dari setiap sel grid

Bo = flowfluid.bO(info.pressure);  % Oil FVF pada tekanan awal
Bw = flowfluid.bW(info.pressure);  % Water FVF pada tekanan awal
Bg = flowfluid.bG(info.pressure);  % Gas FVF pada tekanan awal

% Faktor konversi (mengubah m^3 ke STB jika diperlukan)

conversionFactor_oil = 6.28981;   % m³ ke STB
conversionFactor_gas = 35.3147;   % m³ ke SCF
conversionFactor_water = 6.28981; % m³ ke STB

OWIP = sum(rock.poro .* state0.s(:, 1) .* V_b ./ Bw) * conversionFactor_water;
OOIP = sum(rock.poro .* state0.s(:, 2) .* V_b ./ Bo) * conversionFactor_oil;
OGIP = sum(rock.poro .* state0.s(:, 3) .* V_b ./ Bg) * conversionFactor_gas;

% Menampilkan hasil OOIP
fprintf('Original Water In Place (OWIP): %.2f STB\n', OWIP);
fprintf('Original Oil In Place (OOIP): %.2f STB\n', OOIP);
fprintf('Original Gas In Place (OGIP): %.2f SCF\n', OGIP);
 
pressures = linspace(info.pressure, 100 * psia, 50);

numPressures = length(pressures);
vaporFractions = zeros(1, numPressures);
dew_pressures = zeros(1, numPressures);
dew_vapor = zeros(1, numPressures);

dewCount = 0;
bubbleFound = false;
bubbleIndex = NaN;

for i = 1:numPressures
    p = pressures(i);
    
    % Estimasi K-values (Wilson)
    K = estimateEquilibriumWilson(model.EOSModel, p, info.temp);
    
    % Phase Stability Test
    [stable, x, y] = phaseStabilityTest(model.EOSModel, info.initial, p, info.temp, K);
    
    % Flash Calculation (Rachford-Rice)
    Lfrac = solveRachfordRiceVLE([], K, info.initial); % Liquid fraction
    V = 1 - Lfrac;  % Vapor fraction
    
    vaporFractions(i) = V;
    
    if ~bubbleFound && V > 0
        bubble_pressures = p;
        bubble_vapor = V;
        bubbleIndex = i;
        bubbleFound = true;
    end
    
    if abs(V - 1) < 0.001  % Dew point detection
        dewCount = dewCount + 1;
        dew_pressures(dewCount) = p;
        dew_vapor(dewCount) = V;
    end
    %fprintf('Tekanan: %.2f bar | Vapor Fraction: %.6f | Stable: %d\n', p / barsa, V, stable);
end

dew_pressures = dew_pressures(1:dewCount);
dew_vapor = dew_vapor(1:dewCount);

% Plot Phase Envelope (Pressure vs Vapor Fraction)
figure; hold on;

if ~isnan(bubbleIndex)
    plot(pressures(1:bubbleIndex) / barsa, vaporFractions(1:bubbleIndex), ...
        '-b', 'LineWidth', 1.5, 'DisplayName', 'Liquid Phase Curve');
end

if ~isnan(bubbleIndex) && bubbleIndex < numPressures
    plot(pressures(bubbleIndex:end) / barsa, vaporFractions(bubbleIndex:end), ...
        '-r', 'LineWidth', 1.5, 'DisplayName', 'Two Phase Curve');
end

if ~isnan(bubbleIndex)
    plot(bubble_pressures / barsa, bubble_vapor, 'ro', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Bubble Point');
end

plot(dew_pressures / barsa, dew_vapor, '-x', ...
    'LineWidth', 1.5, 'Color', [0 0.7 0], 'DisplayName', 'Dew Point Curve');

xlabel('Pressure (bar)');
ylabel('Vapor Fraction');
title('Phase Envelope (Pressure vs Vapor Fraction)');
legend('Location', 'best');
grid on;
set(gca, 'YLim', [0 1]);

% Define Wells
W = []; % Inisialisasi array untuk menyimpan sumur

% Sumur Injeksi: INJ1
W = verticalWell(W, G, rock, ...
    1, 1, 10, ...             % Lokasi: baris (I=1), kolom (J=1), perforasi penuh (K=[10])
    'name', 'Inj_1', ...        % Nama sumur
    'compi', [0 0 1], ...       % Produksi hanya minyak
    'val', injRate, ...         % Tekanan dasar sumur (BHP)
    'type', 'rate',...          % Kontrol berdasarkan tekanan dasar
    'sign', 1);                 % Production (Sign = -1) or Injection (Sign = 1)

% Sumur Injeksi: INJ2
W = verticalWell(W, G, rock, ...
    30, 1, 10, ...              % Lokasi: baris (I=30), kolom (J=1), perforasi penuh (K=[10])
    'name', 'Inj_2', ...        % Nama sumur
    'compi', [0 0 1], ...       % Produksi hanya minyak
    'val', injRate, ...         % Tekanan dasar sumur (BHP)
    'type', 'rate',...          % Kontrol berdasarkan tekanan dasar
    'sign', 1);                 % Production (Sign = -1) or Injection (Sign = 1)

% Sumur Injeksi: INJ3
W = verticalWell(W, G, rock, ...
    30, 20, 10, ...             % Lokasi: baris (I=30), kolom (J=20), perforasi penuh (K=[10])
    'name', 'Inj_3', ...        % Nama sumur
    'compi', [0 0 1], ...       % Produksi hanya minyak
    'val', injRate, ...         % Tekanan dasar sumur (BHP)
    'type', 'rate',...          % Kontrol berdasarkan tekanan dasar
    'sign', 1);                 % Production (Sign = -1) or Injection (Sign = 1)

% Sumur Injeksi: INJ4
W = verticalWell(W, G, rock, ...
    1, 20, 10, ...              % Lokasi: baris (I=1), kolom (J=20), perforasi penuh (K=[10])
    'name', 'Inj_4', ...        % Nama sumur
    'compi', [0 0 1], ...       % Produksi hanya minyak
    'val', injRate, ...         % Tekanan dasar sumur (BHP)
    'type', 'rate',...          % Kontrol berdasarkan tekanan dasar
    'sign', 1);                 % Production (Sign = -1) or Injection (Sign = 1)

% Sumur Injeksi: PROD1
W = verticalWell(W, G, rock, ...
    15, 10, 5:10, ...             % Lokasi: baris (I=15), kolom (J=10), perforasi penuh (K=[5:10])
    'name', 'Prod', ...         % Nama sumur
    'compi', [0 1 0], ...       % Produksi hanya minyak
    'val', 100 * barsa, ...     % Tekanan dasar sumur (BHP)
    'type', 'bhp',...           % Kontrol berdasarkan tekanan dasar
    'sign', -1);                % Production (Sign = -1) or Injection (Sign = 1)

for i = 1:numel(W)
    W(i).components = info.injection;
end


% Visualisasi Grid Reservoir
% Membuat struktur data untuk digunakan dalam plotToolbar
data = struct();
data.Porosity = rock.poro;                   % Porositas
data.Kx = rock.perm(:, 1) / (milli * darcy); % Permeabilitas horizontal (Kx)
data.Ky = rock.perm(:, 2) / (milli * darcy); % Permeabilitas horizontal (Ky)
data.Kz = rock.perm(:, 3) / (milli * darcy); % Permeabilitas vertikal (Kz)
data.Pressure = state0.pressure / barsa;     % Tekanan dalam bar
data.SaturationWater = state0.s (:, 1);      % Saturasi air (sw)
data.SaturationOil = state0.s(:, 2);         % Saturasi minyak (sO)
data.SaturationGas = state0.s(:, 3);         % Saturasi gas (sG)

% Gunakan plotToolbar untuk membuat GUI interaktif
figure;
plotToolbar(G, data,'EdgeColor', 'k', 'LineWidth', 0.2); % Menampilkan data reservoir secara interaktif

% Tambahkan pengaturan tampilan 3D
view(2); % Tampilan 3D
axis equal tight;
title('Interactive Reservoir Properties Visualization with Wells');

% Tambahkan colorbar
Colorbar = colorbar;

% Tambahkan sumur ke dalam plot
hold on;
plotWell(G, W, 'color', 'r', 'LineWidth', 2);
hold off;
%% ------------------------------------------------------------------------------------------------------------------------ %%

%% === Phase 4 & 5: Simulation Schedule and Simulate Base Case === %%
% Definisikan timestep untuk produksi dan injeksi menggunakan fungsi ramp-up
dt1 = rampupTimesteps(simTimeProd, nstep, refine); % Timestep untuk periode produksi tanpa injeksi
dt2 = rampupTimesteps(simTimeInj, nstep, refine); % Timestep untuk periode produksi dengan injeksi

% Inisialisasi struktur schedule untuk simulasi
schedule = struct();           % Inisialisasi schedule sebagai struktur kosong
schedule.control = struct([]); % Inisialisasi array control kosong untuk sumur
schedule.step = struct();      % Inisialisasi array step kosong untuk menyimpan timestep

% Gabungkan timestep produksi dan injeksi dalam satu array
schedule.step.val = [dt1; dt2]; % Menyatukan timestep dari dua periode ke dalam schedule

% Duplikasi daftar sumur untuk dua skenario yang berbeda
W1 = W; % Sumur produksi aktif, sumur injeksi nonaktif
W2 = W; % Sumur produksi aktif, sumur injeksi aktif

% Menonaktifkan sumur injeksi hanya pada periode produksi tanpa injeksi
for i = 1:numel(W1)
    if startsWith(W1(i).name, 'Inj')
        W1(i).val = 0;       
        W1(i).type = 'rate';
    end
end

% Sumur injeksi aktif kembali di periode kedua
for i = 1:numel(W2)
    if startsWith(W2(i).name, 'Inj') 
        W2(i).val = injRate;       
        W2(i).type = 'rate';
        W2(i).components = info.injection;
    end
end

% Buat kontrol untuk setiap tahap produksi dan injeksi
schedule.control = repmat(struct('W', []), 2, 1); % Buat array struct dengan field 'W' untuk kontrol
schedule.control(1).W = W1;                       % Kontrol pertama: hanya sumur produksi yang aktif
schedule.control(2).W = W2;                       % Kontrol kedua: sumur produksi dan injeksi aktif

% Mapping control ke step agar setiap timestep menggunakan kontrol yang sesuai
schedule.step.control = [ones(numel(dt1), 1); 2 * ones(numel(dt2), 1)]; % Bagian pertama menggunakan kontrol 1 (hanya produksi), bagian kedua menggunakan kontrol 2 (produksi + injeksi)


nls = NonLinearSolver('useRelaxation', true);

fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true);

[~, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                         'nonlinearsolver', nls, ...
                                         'afterStepFn', fn);

% Inisialisasi array untuk menyimpan viskositas tiap grid dan timestep
numSteps = numel(states);             % Jumlah timestep
numCells = numel(states{1}.pressure); % Jumlah grid cells

viscosityOilGrid = zeros(numCells, numSteps); % Matrix untuk viskositas minyak di setiap timestep
viscosityGasGrid = zeros(numCells, numSteps); % Matrix untuk viskositas gas di setiap timestep
densityOilGrid = zeros(numCells, numSteps);
densityGasGrid = zeros(numCells, numSteps);


% Iterasi pada setiap timestep untuk menghitung viskositas di setiap grid
for t = 1:numSteps
    state = states{t};
    
    % Ambil data dari states
    pressure = state.pressure;
    temperature = state.T;
    x = state.x;
    y = state.y;
    Z_L = state.Z_L;
    Z_V = state.Z_V;
    rho = state.rho;

    % Konversi x dan y ke cell array untuk kompatibilitas dengan computeViscosity
    x_cell = num2cell(x, 2);
    y_cell = num2cell(y, 2);

    % Inisialisasi viskositas per grid
    muL = zeros(numCells, 1);
    muV = zeros(numCells, 1);
    rhoOil = zeros(numCells, 1);
    rhoGas = zeros(numCells, 1);

    % Iterasi setiap sel grid untuk hitung viskositas
    for cell = 1:numCells
        P = pressure(cell);
        T = temperature(cell);
        Zl = Z_L(cell);
        Zv = Z_V(cell);
        compL = x_cell{cell};
        compV = y_cell{cell};

        % Hitung viskositas minyak
        if sum(compL) > 0
            muL(cell) = computeViscosity(model.EOSModel.PropertyModel, model.EOSModel, P, compL, Zl, T, true);
        else
            muL(cell) = NaN;
        end

        % Hitung viskositas gas
        if sum(compV) > 0
            muV(cell) = computeViscosity(model.EOSModel.PropertyModel, model.EOSModel, P, compV, Zv, T, false);
        else
            muV(cell) = NaN;
        end

        rhoOil(cell) = rho(cell, 2);
        rhoGas(cell) = rho(cell, 3);
    end

    % Simpan viskositas untuk setiap grid pada timestep ini
    viscosityOilGrid(:, t) = muL;
    viscosityGasGrid(:, t) = muV;
    densityOilGrid(:, t) = rhoOil;
    densityGasGrid(:, t) = rhoGas;
end

figure;
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(2);
axis equal tight;
grid on;

for t = 1:numSteps
    cla;
    muL_grid = viscosityOilGrid(:, t);
    plotCellData(G, muL_grid);
    title(sprintf('Distribusi Viskositas Minyak, Timestep %d', t));
    pause(0.1);
end

figure;
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(2);
axis equal tight;
grid on;

for t = 1:numSteps
    cla;
    muV_grid = viscosityGasGrid(:, t);
    plotCellData(G, muV_grid);
    title(sprintf('Distribusi Viskositas Gas, Timestep %d', t));
    pause(0.1);
end

figure;
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(2);
axis equal tight;
grid on;

for t = 1:numSteps
    cla;
    densityOil_grid = densityOilGrid(:, t);
    plotCellData(G, densityOil_grid);
    title(sprintf('Distribusi Densitas Oil, Timestep %d', t));
    pause(0.1);
end

figure;
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(2);
axis equal tight;
grid on;

for t = 1:numSteps
    cla;
    densityGas_grid = densityGasGrid(:, t);
    plotCellData(G, densityGas_grid);
    title(sprintf('Distribusi Densitas Gas, Timestep %d', t));
    pause(0.1);
end

% Ekstrak data wellSols dari setiap langkah waktu dalam states
wellSols = cellfun(@(x) x.wellSol, states, 'UniformOutput', false);

% Daftar properti yang ingin diekstrak
properties = {'status', 'type', 'bhp', 'qWs', 'qOs', 'qGs', 'qSs', ...
              'qs', 'qTs', 'qTr', 'qGr', 'qOr', 'qWr', 'wcut', 'gcut', ...
              'ocut', 'gor'};

% Inisialisasi struktur untuk menyimpan hasil
numProdWells = 1; % Jumlah sumur produksi
numInjWells = 4;  % Jumlah sumur injeksi

wellDataProd = struct();
wellDataInj = struct();

% Ekstrak tanda sumur (sign) dari setiap langkah waktu
signs = cellfun(@(step) arrayfun(@(well) well.sign, step, 'UniformOutput', true), wellSols, 'UniformOutput', false);

% Loop untuk mengekstrak semua data dari wellSols dan memisahkan berdasarkan tanda
for t = 1:numel(wellSols)               % Iterasi untuk setiap langkah waktu
    % Filter sumur produksi dan injeksi
    prodIndices = find(signs{t} == -1); % Indeks untuk sumur produksi
    injIndices = find(signs{t} == 1);   % Indeks untuk sumur injeksi

    % Pastikan ada cukup sumur dalam data
    if length(prodIndices) < numProdWells
        error('Jumlah sumur produksi yang ditemukan kurang dari yang diharapkan.');
    end
    if length(injIndices) < numInjWells
        error('Jumlah sumur injeksi yang ditemukan kurang dari yang diharapkan.');
    end

    % Simpan data untuk setiap sumur produksi
    for w = 1:numProdWells
        for p = 1:numel(properties)
            propName = properties{p};
            if t == 1
                wellDataProd(w).(propName) = {};
            end
            wellDataProd(w).(propName){t,1} = wellSols{t}(prodIndices(w)).(propName);
        end
    end

    % Simpan data untuk setiap sumur injeksi
    for w = 1:numInjWells
        for p = 1:numel(properties)
            propName = properties{p};
            if t == 1
                wellDataInj(w).(propName) = {};
            end
            wellDataInj(w).(propName){t,1} = wellSols{t}(injIndices(w)).(propName);
        end
    end
end

%% Bagian Konversi Unit
% Faktor konversi dari SI ke field unit
SI_to_STB_per_day = 86400 / 0.159;    % Konversi m³/s ke STB/d (untuk minyak & air)
SI_to_SCF_per_day = 86400 / 28.3168;  % Konversi m³/s ke SCF/d (untuk gas)
SI_to_psi = 0.000145038;              % Konversi Pascal (Pa) ke psi

% Properti yang perlu dikonversi
convertProps = {'qWs', 'qOs', 'qSs', 'qs', 'qTs', 'qTr', 'qOr', 'qWr'};
convertPropsGas = {'qGs', 'qGr'};
convertPropsPressure = {'bhp'};  % Properti tekanan yang perlu dikonversi ke psi

% Struktur baru untuk menyimpan data dalam field unit
wellDataProd_Field = wellDataProd;
wellDataInj_Field = wellDataInj;

% Lakukan konversi hanya untuk properti yang membutuhkan
for w = 1:numProdWells
    for p = 1:numel(properties)
        propName = properties{p};
        if ismember(propName, convertProps)
            wellDataProd_Field(w).(propName) = cellfun(@(x) abs(x) * SI_to_STB_per_day, wellDataProd(w).(propName), 'UniformOutput', false);
        elseif ismember(propName, convertPropsGas)
            wellDataProd_Field(w).(propName) = cellfun(@(x) abs(x) * SI_to_SCF_per_day, wellDataProd(w).(propName), 'UniformOutput', false);
        elseif ismember(propName, convertPropsPressure)
            wellDataProd_Field(w).(propName) = cellfun(@(x) x * SI_to_psi, wellDataProd(w).(propName), 'UniformOutput', false);
        end
    end
end

for w = 1:numInjWells
    for p = 1:numel(properties)
        propName = properties{p};
        if ismember(propName, convertProps)
            wellDataInj_Field(w).(propName) = cellfun(@(x) x * SI_to_STB_per_day, wellDataInj(w).(propName), 'UniformOutput', false);
        elseif ismember(propName, convertPropsGas)
            wellDataInj_Field(w).(propName) = cellfun(@(x) x * SI_to_SCF_per_day, wellDataInj(w).(propName), 'UniformOutput', false);
        elseif ismember(propName, convertPropsPressure)
            wellDataInj_Field(w).(propName) = cellfun(@(x) x * SI_to_psi, wellDataInj(w).(propName), 'UniformOutput', false);
        end
    end
end

% Nama file Excel untuk menyimpan data
outputFilename = 'WellData-FieldCase3-1.xlsx';

% Menyimpan data produksi ke dalam Excel (setiap sumur pada sheet yang berbeda)
for w = 1:numProdWells
    prodTable = struct2table(wellDataProd_Field(w));
    writetable(prodTable, outputFilename, 'Sheet', ['Prod' num2str(w)], 'WriteVariableNames', true);
end

% Menyimpan data injeksi ke dalam Excel (setiap sumur pada sheet yang berbeda)
for w = 1:numInjWells
    injTable = struct2table(wellDataInj_Field(w));
    writetable(injTable, outputFilename, 'Sheet', ['Inj' num2str(w)], 'WriteVariableNames', true);
end
%% ------------------------------------------------------------------------------------------------------------------------ %%

%% === Phase 6: Visualization and Plot Reservoir Data === %%

% Menampilkan hasil simulasi terakhir
i = numel(states); % Ambil langkah waktu terakhir

figure;
% Plot 1: Tekanan (Pressure)
subplot(3, 1, 1);
plotCellData(G, states{i}.pressure, 'EdgeColor', 'none');
axis tight; set(gca, 'xticklabel', [], 'yticklabel', []);
colorbar;
title('Tekanan Reservoir');
ylabel('Pressure [bar]');
view(90, 0);

% Plot 2: Saturasi Gas (Gas Saturation)
subplot(3, 1, 2);
plotCellData(G, states{i}.s(:, 3), 'EdgeColor', 'none'); % Kolom 3 adalah gas
axis tight; set(gca, 'xticklabel', [], 'yticklabel', []);
clim([0, 1]); colorbar;
title('Saturasi Gas');
ylabel('Gas Saturation');
view(90, 0);

% Plot 3: Fraksi Mole CO₂ (CO2 Mole Fraction)
subplot(3, 1, 3);
plotCellData(G, states{i}.components(:, 2), 'EdgeColor', 'none'); % Komponen CO₂
axis tight; set(gca, 'xticklabel', [], 'yticklabel', []);
clim([0, 1]); colorbar;
title('Fraksi Mole CO₂');
ylabel('CO₂ Mole Fraction');
view(90, 0);