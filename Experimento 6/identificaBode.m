posMax = []; %vetor para posi��o m�xima
posMin = []; %vetor para posi��o m�nima

listaOmega = [25 26 27 28 29 30 31 32 33];

for omega=listaOmega
    
    fprintf('\nIniciando Simulacao para Omega = %i\n', omega);
	% vetor para posi��o 
	pos = [];
   
    % inicializa a remoteAPI
    vrep=remApi('remoteApi');
    % por seguranca, fecha todas as conexoes abertas
    vrep.simxFinish(-1);
    %se conecta ao V-rep usando a porta 19997 (mesmo computador)
    clientID=vrep.simxStart('127.0.0.1',19997,true,true,5000,5);

    %d� play na simula��o
    vrep.simxStartSimulation(clientID,vrep.simx_opmode_blocking);


    if clientID == -1
        fprintf('aconteceu algum problema na conexao com o servidor da remoteAPI!\n');
        return;
    end

    % escolhe a comunicacao sincrona
    vrep.simxSynchronous(clientID,true);

    % pega um handle para a massa superior (chamada de topPlate). Ele sera
    % usado para aplicarmos a forca.
    [err,h]=vrep.simxGetObjectHandle(clientID,'SpringDamper_topPlate',vrep.simx_opmode_oneshot_wait);

    % pega a posicao inicial da massa
    [res,posTopPlate]=vrep.simxGetObjectPosition(clientID,h,-1,vrep.simx_opmode_oneshot_wait);

    % coloca a posicao (tres coordenadas, mas apenas a z vai interessar) na
    % variavel position.
    position = posTopPlate;

    % espera sincronizar
    vrep.simxSynchronousTrigger(clientID);

    % ajusta o passo de tempo (deve ser o mesmo que esta ajustado no v-rep)
    timeStep = 0.01;

    %inicializa o tempo da simulacao.
    timeSim=0;
    % laco principal (fica ate nao haver mudanca significativa na posicao da
    % massa)

    fprintf('posicao inicial = %.7f\n',position(1,end));
    while timeSim(1,end) < 3
        % aplica a forca senoidal
        [res retInts retFloats retStrings retBuffer]=vrep.simxCallScriptFunction(clientID,'myFunctions',vrep.sim_scripttype_childscript,'addForceTo',[h],[0,0,0,0,0,-1*sin(omega*timeSim(1,end))],[],[],vrep.simx_opmode_blocking);

        % espera sincronizar
        vrep.simxSynchronousTrigger(clientID);
        timeSim=[timeSim timeSim(1,end)+timeStep];

        % faz a leitura da posicao da massa
        [res,posTopPlate]=vrep.simxGetObjectPosition(clientID,h,-1,vrep.simx_opmode_oneshot_wait);

        %guarde a nova posicao posTopPlate aqui...
        pos(end+1) = posTopPlate(3);

        %imprime na tela
        fprintf('t=%.3f -> pos = [%.5f,%.5f,%.5f]\n',timeSim(1,end),posTopPlate(1),posTopPlate(2),posTopPlate(3));

    end

    posMax(end+1) = max(pos); % adiciona posi��o m�xima no vetor posMax
    posMin(end+1) = min(pos); % adiciona posi��o m�nima no vetor posMin

    % chama o destrutor
    vrep.delete(); % call the destructor!

    fprintf('Simulacao terminada!\n');
    %d� stop na simula��o
    vrep.simxStopSimulation(clientID,vrep.simx_opmode_blocking);

end

%Amplitude
Mp = (posMax - posMin)/(2*0.0049);

%plot de Mp em fun��o da frequ�ncia
plot(listaOmega, Mp,'b',listaOmega, Mp,'r o');
hold on;
title("M�ximo Sobressinal em fun��o da frequ�ncia");
xlabel("Omega");ylabel("Mp");
set(gca,'xscale','log');
grid();
    
%plot de Mp dB em fun��o da frequ�ncia
Mp_dB = mag2db(Mp);
plot(listaOmega, Mp_dB, 'b', listaOmega, Mp_dB, 'r o');
hold on;
title("Diagrama de Bode para o M�ximo Sobressinal");
xlabel("Omega");ylabel("Mp [dB]");
set(gca,'xscale','log')
grid();


%Plot da fun��o de transfer�ncia

%Par�metros
k = 177.79;
b = 2.053;
m = 0.16388;
G = tf([0.0049*k], [m b k]);


figure; Mp_norm = Mp_dB-46.2;
semilogx(listaOmega, Mp_norm, 'ro');legend('Experimental');
hold on;
bodeplot(G);legend('te�ricos')
grid;
xlabel('\omega');ylabel('G(j\omega)');
title('Diagrama de Bode: Valores experimentais e te�ricos')






