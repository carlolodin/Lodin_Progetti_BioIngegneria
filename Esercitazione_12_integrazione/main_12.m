%% Integrazione
x = 0:0.01:10;
f = @(x) 7 + sin(12 * pi * x);
y = f(x);

area_trapz = trapz(x, y);
area_rect_pm = 0;
area_rect_l = 0;
area_rect_r = 0;

for i = 0:length(x)-1
    area_rect_pm = area_rect_pm + f(i*0.01 + 0.005) * 0.01; % punto medio
    area_rect_l = area_rect_l + f(i*0.01 + 0.005) * 0.01; %punto sinistro
    area_rect_r = area_rect_r + f((i+1)*0.01 + 0.005) * 0.01; % punto destro
end

disp(['Area calcolata con trapz: ', num2str(area_trapz)]);
disp(['Area calcolata con rettangoli (punto medio): ', num2str(area_rect_pm)]);
disp(['Area calcolata con rettangoli (punto sinistro): ', num2str(area_rect_l)]);
disp(['Area calcolata con rettangoli (punto destro): ', num2str(area_rect_r)]);

%% Derivazione
f_prime = @(x) 12 * pi * cos(12 * pi * x);
f_prime_forward = zeros(size(x)-1);
f_prime_backward = zeros(size(x)-1);
f_prime_central = zeros(size(x)-2);

for i = 1:length(x)-1
    f_prime_forward(i) = (y(i+1) - y(i)) / 0.01;
    f_prime_backward(i) = (y(i+1) - y(i)) / 0.01;
    if i < length(x)-1 && i > 1
        f_prime_central(i) = (y(i+1) - y(i-1)) / 0.02;
    end
end

figure;
hold on;
plot(x, f_prime(x), 'k', 'DisplayName', 'Derivata Analitica');
plot(x(1:end-1), f_prime_forward, 'r', 'DisplayName', 'Forward Difference');
plot(x(2:end), f_prime_backward, 'g', 'DisplayName', 'Backward Difference');
plot(x(2:end-1), f_prime_central, 'b', 'DisplayName', 'Central Difference');
xlabel('x');
ylabel('f''(x)');
title('Derivate di f(x)');
legend show;
xlim([0, 1]);

%% Metodo di eulero esplicito e runge kutta
dt = 0.01;
y_euler = zeros(size(x));
y_euler(1) = y(1);  % Condizione iniziale
for i = 1:length(x)-1
    y_euler(i+1) = y_euler(i) + dt * f_prime(i * dt);
end

y_rk = zeros(size(x));
y_rk(1) = y(1);  % Condizione iniziale
for i = 1:length(x)-1
    k1 = f_prime(i * dt);
    k2 = f_prime(i * dt + 0.5 * dt);
    k3 = f_prime(i * dt + 0.5 * dt);
    k4 = f_prime((i + 1) * dt);
    y_rk(i+1) = y_rk(i) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

figure;
hold on;
plot(x, f(x), 'k', 'DisplayName', 'f(x)');
plot(x, y_euler, 'r', 'DisplayName', 'Euler Esplicito');
plot(x, y_rk, 'g', 'DisplayName', 'Runge Kutta');
xlabel('x');
ylabel('y');
title('Soluzioni di f(x) con Euler e Runge Kutta');
legend show;
xlim([0, 3]);


