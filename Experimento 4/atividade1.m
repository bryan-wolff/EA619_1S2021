% Atividade 1

close all; clear all; clc;

% Equação diferencial: 
%     y'''(t) + (2+6p)y''(t) + (9+12p)y'(t) + 18y(t) = 18x(t)
% em que p é um parâmetro desconhecido, pertencente à faixa p ? [0.1,1.2]

%---------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------%

% Questão 3
p = linspace(0.1,1.2,50);

% Script para a chamada do diagrama desenvolvido no item 1
for i = 1:50
    simOut = sim('sistemaLinearOrdem3','SrcWorkspace','current','maxstep','0.1');
end

%---------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------%

%Questão 4

% (a)
fprintf('Item (a)\n');
% valor final em que y estabiliza
y_f = simOut.yout{1}.Values.Data(end);

% [maior valor de y, posição do maior valor de y]
[m,pm] = max(simOut.yout{1}.Values.Data);

% [valor máximo para o sobressinal, índice do valor máximo]
[y_max, i] = max(m);
m_max = (y_max - y_f)/y_f;
fprintf("Valor máximo para o sobressinal (Mp): %.3f \n", m_max);

% valor de p para Mp máximo
p_max = p(i);
fprintf('Valor de p para o máximo Mp: %.3f \n', p_max);

% valor de pm para Mp máximo
pm_max = pm(i);
t_mp_max = simOut.tout(pm_max);
fprintf('Instante de tempo para Mp máximo: %.3f s\n', t_mp_max);

% t_s: tempo de subida
% valor de referência na entrada: degrau unitário (= 1)

len_simOut = length(simOut.yout{1}.Values.Data);
for k = 1:len_simOut
    y(k) = simOut.yout{1}.Values.Data(k,i);
end

for j = 1:len_simOut
    j_index = j;
    if y(j_index) > 1
        break
    end
end
t_s = simOut.tout(j_index);
fprintf('Tempo de subida: %.3f s \n', t_s);

% t_a: tempo de acomodação
% valor de referência na entrada: degrau unitário (= 1)
% (1 - 0.02*1) = 0.98
% (1 + 0.02*1) = 1.02

k_in = 0;
k_aux = 0;
for k = 1:len_simOut
    if (y(k) >= (1-0.02*1)) && (y(k) <= (1+0.02*1)) 
        if (k_aux == 0)
            k_in = k;
            k_aux = 1;
        end
    else
        k_aux = 0;
    end
end
t_a = simOut.tout(k_in);
fprintf('Tempo de acomodação: %.3f s\n', t_a);

fprintf('------------------------------------------------------\n');

%---------------------------------------------------------------------------------%
% (b)
fprintf('Item (b)\n');
% mp < 0.07 --> mp = (m - y_f)/y_f < 0.07
contador = 1;
p7 = [];
for t = 1:50
  if (((m(t) - y_f)/y_f) <= 0.07)
      p7(contador) = p(t);
      contador = contador + 1;
  end
end

fprintf('Subfaixa de p para Mp <= 0.07: (%.3f , %.3f)\n', p7(1),p7(end));

fprintf('------------------------------------------------------\n');

%---------------------------------------------------------------------------------%
% (c)
fprintf('Item (c)\n');

% função para definir o p que gera Mp mais próximo de 7%
prox = 1;
for n = 1:50
  if (abs(((m(n) - y_f)/y_f) - 0.07) < prox) %se diferença for menor que a anterior
      prox = abs(((m(n) - y_f)/y_f) - 0.07); %redefine prox
      p_prox = p(n); %redefine p que gera o sinal
      mp_p_prox = m(n) - y_f;
  end
end

% função de transferência
func_transf = tf(18,[1 (2+6*p_prox) (9+12*p_prox) 18]);
polos = pole(func_transf);
fprintf('Polos da função: \n');
disp(polos);

fprintf('------------------------------------------------------\n');

%---------------------------------------------------------------------------------%
% (d)
fprintf('Item (d)\n');

contador = 1;
p1 = [];
for t = 1:50
  if (((m(t) - y_f)/y_f) <= 0.01)
      p1(contador) = p(t);
      contador = contador + 1;
  end
end

fprintf('Subfaixa de p para Mp <= 0.01: (%.3f , %.3f)\n', p1(1),p1(end));



