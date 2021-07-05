% Atividade 2
clear all; close all; clc;

% p^2theta + g/l*sen(theta) = 0
% x1 = theta
% x2 = theta'
% A = [df1/x1  df1/x2] = [0              1]   B = [0]   u = 0
%     [df2/x1  df2/x2]   [(-g/l)cos(x1)  0]       [0]

% vari�veis globais
global g; % acelera��o da gravidade
global l; % comprimento da haste

l = 1;
g = 9.81;

tempoTotal = 10;
optOde = odeset('maxStep', 0.05);

% condi��es iniciais
x0 = [pi/10 0];

% sistema linear
outOde_linear10 = ode45(@linSys, [0 tempoTotal], x0, optOde);

% sistema n�o linear
outOde_nao_linear10 = ode45(@notlinSys, [0 tempoTotal], x0, optOde);

% plot das fun��es
figure();
plot(outOde_linear10.x, outOde_linear10.y);
hold on
plot(outOde_nao_linear10.x, outOde_nao_linear10.y)
xlabel t(s); grid;
legend("\theta_{l}","\omega_{l}","\theta_{nl}","\omega_{nl}")
title("\theta e \omega = \theta' em fun��o do tempo para condi��o inicial (\pi/10, 0)")

x02 = [pi/4 0];
% sistema linear
outOde_linear4 = ode45(@linSys, [0 tempoTotal], x02, optOde);

% sistema n�o linear
outOde_nao_linear4 = ode45(@notlinSys, [0 tempoTotal], x02, optOde);

% plot das fun��es
figure();
plot(outOde_linear4.x, outOde_linear4.y);
hold on
plot(outOde_nao_linear4.x, outOde_nao_linear4.y)
xlabel t(s); grid;
legend("\theta_{l}","\omega_{l}","\theta_{nl}","\omega_{nl}")
title("\theta e \omega = \theta' em fun��o do tempo para condi��o inicial (\pi/4,0)")

% plano de fases para a primeira condi��o inicial
figure()
plot(outOde_nao_linear10.y(1,:),outOde_nao_linear10.y(2,:),'r');
title('Plano de fase para condi��o inicial (\pi/10,0)');
hold on
plot(outOde_linear10.y(1,:),outOde_linear10.y(2,:),'b');
xlabel('\theta (rad)');
ylabel('\omega (rad/s)'); grid;
legend('N�o linear','Linear');

% plano de fases para a segunda condi��o inicial
figure()
plot(outOde_nao_linear4.y(1,:),outOde_nao_linear4.y(2,:),'r');
title('Plano de fase para condi��o inicial (\pi/4,0)');
hold on
plot(outOde_linear4.y(1,:),outOde_linear4.y(2,:),'b');
xlabel('\theta (rad)');
ylabel('\omega (rad/s)'); grid;
legend('N�o linear','Linear');

% gerando o video 1
% geraVideo(outOde_linear10.x,outOde_nao_linear10.y,outOde_linear10.y,tempoTotal,l,'pendSimplesCondicao1');

% gerando o video 2
% geraVideo(outOde_linear4.x,outOde_nao_linear4.y,outOde_linear4.y,tempoTotal,l,'pendSimplesCondicao2');

% Defasagem
defasagem10 = dif(outOde_nao_linear10, outOde_linear10);
title('Defasagem entre os modelos linear e n�o linear para \theta_0 = \pi/10');
defasagem4 = dif(outOde_nao_linear4, outOde_linear4);
title('Defasagem entre os modelos linear e n�o linear para \theta_0 = \pi/4');

disp("Intervalos de tempo com defasagem superior a 2� para pi/10")
grau_2(defasagem10,outOde_linear10.x)
fprintf("\n\n")
disp("Intervalos de tempo com defasagem superior a 2� para pi/4")
grau_2(defasagem4,outOde_linear4.x)
fprintf("\n")

function defasagem = dif(nao_linear,linear)
    defasagem = abs(rad2deg(nao_linear.y(1,:) - linear.y(1,:)));
    figure()
    plot(nao_linear.x, defasagem); grid;
    ylabel('Defasagem (�)'); xlabel('Tempo (s)');
end

function intervalo_tempo = grau_2(defasagem, nao_linear)
    %Encontra os intervalos temporais com defasagem superiores a 2 graus
    aux = 0;
    for j = 1:length(defasagem)
        if (defasagem(1,j) >= 2) && (aux == 0)
            aux = 1;
            fprintf("T: [%.3f ,", nao_linear(j))
        elseif (defasagem(1,j) <= 2) && (aux == 1)
            aux = 0;
            fprintf(" %.3f] s\n", nao_linear(j))
        end
    end
end

function dX = linSys(t,x)
global g;
global l;
g = 9.81;
l = 1;

A = [0 1; (-g/l)*cos(x(1)) 0];
B = [0;0];
u = 0;

dX = A*x + B*u;
end

function dX_nao_linear = notlinSys(t,x)
global g;
global l;
g = 9.81;
l = 1;

dX_nao_linear(1,1) = x(2);
dX_nao_linear(2,1) = -(g/l)*sin(x(1));
end

function geraVideo(tempo,estadosN,estadosL,tempoTotal,L,titulo)
% entradas:
% tempo: vetor contendo os instantes de tempo da simula ?c~ao
% estadosN: vetor contando os estados ( theta e dot_theta )
% da simula��o n�o linear
% estadosL: vetor contando os estados ( theta e dot_theta )
% da simula��o linear
% tempoTotal: tempo total da simula��o
% L: comprimento do fio ( metros )
writerObj = VideoWriter (titulo,'Motion JPEG AVI');
writerObj.FrameRate = ceil(size(tempo,2) / tempoTotal);
open (writerObj);
fig = figure();
ori = [0 0]; % origem do p�ndulo
f = 1;
while f <= length(tempo)
    theta_nao_linear = estadosN(1,f);
    % determine aqui as coordenadas do corpo para o caso n�o linear (use L)
    x_nao_linear = L*sin(theta_nao_linear);
    y_nao_linear = -L*cos(theta_nao_linear);
    
    % determine aqui as coordenadas do corpo para o caso linearizado(use L)
    theta_linear = estadosL(1,f);
    x_linear = L*sin(theta_linear);
    y_linear = -L*cos(theta_linear);

    axis ([-2 2 -2.5 0.5]);
    hold on ;
    
    % desenhe os p�ndulos aqui
    %n�o-linear
    plot(x_nao_linear,y_nao_linear,'o','MarkerSize',10,'MarkerFaceColor','r')
    line([0 x_nao_linear],[0 y_nao_linear]);
    
    %linear
    plot(x_linear,y_linear,'o','MarkerSize',10,'MarkerFaceColor','b')
    line([0 x_linear],[0 y_linear]);
    
    legend('Massa NL','Haste NL','Massa L','Haste L'); grid;
    
    hold off;
    F = getframe;
    writeVideo(writerObj,F);
    clf(fig);
    f = f + 1;
end
close(writerObj);
end











