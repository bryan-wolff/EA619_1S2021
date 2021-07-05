% Atividade 3
clear all; close all; clc;

global g;
global L1;
global L2;
global m1;
global m2;

% arbitrando valores para as variáveis
g = 9.81; % aceleração da gravidade
L1 = 1; % comprimento da haste 1
L2 = 0.5; % comprimento da haste 2
m1 = 2; % corpo de massa 1
m2 = 0.5; % corpo de massa 2

% condição inicial do sistema
x0 = [pi/4 0 pi/4 0];

tempoTotal = 10;
optOde = odeset('maxStep', 0.05);

% solução a partir da função de pêndulo duplo
outOde = ode45(@penduloDuplo,[0 tempoTotal],x0,optOde);

% plot da função
figure();
plot(outOde.x, outOde.y);
xlabel t(s); grid;
legend("\theta_{1}","\omega_{1}","\theta_{2}","\omega_{2}");
title("Pêndulo Duplo: Comportamento de \theta e \omega em função do tempo");

% planos de fase
figure();
plot(outOde.y(1,:),outOde.y(2,:));
xlabel('\omega_1[rad/s]'); grid;
ylabel('\theta_1[rad]');
title("Pêndulo duplo: Plano de fase para \theta_{1} e \omega_{1}");

figure()
plot(outOde.y(3,:),outOde.y(4,:)); xlabel t;
xlabel('\omega_2[rad/s]'); grid;
ylabel('\theta_2[rad]');
title("Pêndulo duplo: Plano de fase para \theta_{2} e \omega_{2}");

% gerando video
geraVideo(outOde.x,outOde.y,tempoTotal,L1,L2,'penduloDuplo.avi');

function dX = penduloDuplo(t,x)
global L1;
global L2;
global m1;
global m2;
global g;

% variaveis de estado x1 = theta_1 , x2 = w_1 , x3 = theta_2 , x4 = w_2
dX(1,1) = x(2);
num1 = -g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3)));
den1 = L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)));
num2 = 2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3)));
den2 = L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3)));
dX(2,1) = num1/den1;
dX(3,1) = x(4);
dX(4,1) = num2/den2;
end

function geraVideo(tempo,estadosN,tempoTotal,L1,L2,titulo)
% entradas :
% tempo : vetor contendo os instantes de tempo da simula ?c~ao
% estadosN : vetor contando os estados ( theta e dot_theta )
%da simulação não linear
% estadosL : vetor contando os estados ( theta e dot_theta )
%da simulação linear
% tempoTotal : tempo total da simula ?c~ao
% L: comprimento do fio ( metros )
writerObj = VideoWriter(titulo,'Motion JPEG AVI');
writerObj.FrameRate = ceil(size(tempo,2)/tempoTotal);
open(writerObj);
fig = figure();
ori = [0 0]; % origem do p^endulo
f = 1;
while f <= length(tempo)
    theta1 = estadosN(1,f);
    theta2 = estadosN(3,f);
    % determine aqui as coordenadas do corpo para o corpo 1
    x1 = L1*sin(theta1);
    y1 = -L1*cos(theta1);
    % determine aqui as coordenadas do corpo para o corpo 2
    x2 = x1 + L2*sin(theta2);
    y2 = y1 - L2*cos(theta2);

    axis([-2 2 -2.5 0.5]);
    hold on;
    % desenhe os pêndulos aqui
    plot(x1,y1,'o','MarkerSize',10,'MarkerFaceColor','r')
    line([0 x1], [0 y1])
    plot(x2,y2,'o','MarkerSize',10,'MarkerFaceColor','b')
    line([x1 x2], [y1 y2])
    hold off;
    F = getframe;
    writeVideo(writerObj,F);
    clf(fig);
    f = f + 1;
end
close(writerObj);
end



