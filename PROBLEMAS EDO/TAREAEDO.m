% Script de Simulación Numérica: Oscilador de Van der Pol mediante RK4

% --- [I] IDENTIFICAR: Parámetros del modelo y condiciones iniciales ---
mu = 1.5;            % Coeficiente de amortiguamiento no lineal
h = 1.0;            % Nuestro tamaño de paso temporal (s)
t = 0:h:20;          % Vector temporal: simularemos 20 segundos
N = length(t);       % Total de puntos a calcular

x_pos = zeros(1, N); % Memoria para el desplazamiento estructural
v_vel = zeros(1, N); % Memoria para la velocidad

% Arrancamos asumiendo el desplazamiento máximo y reposo instantáneo
x_pos(1) = 2.0;        
v_vel(1) = 0.0;        

% [P] PLANTEAR: Sistema de ecuaciones de primer orden
% Ecuación 1: La derivada de la posición es la velocidad
dx_dt = @(t, x_pos, v_vel) v_vel;

% Ecuación 2: Aceleración (despejada de la ecuación de Van der Pol)
dv_dt = @(t, x_pos, v_vel) mu * (1 - x_pos^2) * v_vel - x_pos;

% [E] EJECUTAR: Integración numérica con algoritmo RK4
for i = 1:(N-1)
    % 1. Evaluamos pendientes al inicio del paso
    k1_x = dx_dt(t(i), x_pos(i), v_vel(i));
    k1_v = dv_dt(t(i), x_pos(i), v_vel(i));

    % 2. Evaluamos en el primer punto medio
    k2_x = dx_dt(t(i) + h/2, x_pos(i) + (h/2)*k1_x, v_vel(i) + (h/2)*k1_v);
    k2_v = dv_dt(t(i) + h/2, x_pos(i) + (h/2)*k1_x, v_vel(i) + (h/2)*k1_v);

    % 3. Evaluamos en el segundo punto medio
    k3_x = dx_dt(t(i) + h/2, x_pos(i) + (h/2)*k2_x, v_vel(i) + (h/2)*k2_v);
    k3_v = dv_dt(t(i) + h/2, x_pos(i) + (h/2)*k2_x, v_vel(i) + (h/2)*k2_v);

    % 4. Evaluamos al final del intervalo
    k4_x = dx_dt(t(i) + h, x_pos(i) + h*k3_x, v_vel(i) + h*k3_v);
    k4_v = dv_dt(t(i) + h, x_pos(i) + h*k3_x, v_vel(i) + h*k3_v);

    % Actualizamos los estados ponderando las 4 pendientes
    x_pos(i+1) = x_pos(i) + (h/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
    v_vel(i+1) = v_vel(i) + (h/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
end

% [E] EVALUAR: Salida de datos en consola y visualización 

% Mostramos la tabla con espaciado fijo para mantener la alineación estricta
fprintf('\n%-12s %-12s %-18s %-15s\n', 'Iteración', 'Tiempo(s)', 'Desplazamiento', 'Velocidad');
fprintf('--------------------------------------------------------------\n');
for i = 1:11 
    fprintf('%-12d %-12.2f %-18.4f %-15.4f\n', i-1, t(i), x_pos(i), v_vel(i));
end

fprintf('\n>>> ESTADO FINAL (t = %.2f s):\n', t(end));
fprintf('Desplazamiento: %.4f\n', x_pos(end));
fprintf('Velocidad: %.4f\n\n', v_vel(end));

% Definición de paleta de colores personalizada (Códigos Hexadecimales)
color_pos = '#0072BD';  % Azul corporativo
color_vel = '#D95319';  % Naranja intenso
color_fase = '#7E2F8E'; % Púrpura oscuro

% Gráfica 1: Evolución temporal
figure(1);
plot(t, x_pos, 'Color', color_pos, 'LineWidth', 2); hold on;
plot(t, v_vel, 'Color', color_vel, 'LineWidth', 2);
grid on;
title('Evolución Temporal de la Estructura (Van der Pol)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Desplazamiento x(t)', 'Velocidad v(t)', 'Location', 'best');

% Gráfica 2: Diagrama de espacio de fases
figure(2);
plot(x_pos, v_vel, 'Color', color_fase, 'LineWidth', 1.5);
grid on;
title('Diagrama de Espacio de Fases (Ciclo Límite)');
xlabel('Desplazamiento x(t)');
ylabel('Velocidad v(t)');