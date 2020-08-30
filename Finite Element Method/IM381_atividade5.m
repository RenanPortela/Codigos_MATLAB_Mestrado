%% Program: Área, Centróide e Momento de inércia
%
%   @DESCRIÇÃO Atividade 5: Definir a área da geometria, seu centróide e
%   momento de inércia
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar memória e fechar janelas
clear all
close all
clc

%% Matriz de Coordenadas do sistema 
%   coord=[nº | X | Y] => Matriz de coordenadas dos nós

nnos = 5; %número de nós menos 1
nel = nnos - 1;

A = 0; %somatório das áreas dos elementos
e_Area = zeros(6,1);
e_Centroide = zeros(6,1);
e_Inercia = zeros(6,1);

while norm(pi/4-A)>1e-4
    A = 0;
    coord = zeros(nnos,3);

    coord(1,1) = 1;
    coord(1,2) = 0;
    coord(1,3) = 0;

    theta = (pi/2)/(nnos-1);

    for i = 1:nnos
    coord(i+1,1) = i+1; %número do nó
    coord(i+1,2) = cos(theta*(i-1))+coord(1,2); %coordenada x
    coord(i+1,3) = sin(theta*(i-1))+coord(1,3); %coordenada y
    end



    %-------Matriz de incidência ---------------------------------------------%
    % inci=[nº | tab.mat | tab.geo | nó1 | nó2 ] => Matriz de
    % incidência do elemento

    inci = zeros(nnos-1,3); 

    for i = 1:nnos-1 
        inci(i,1) = 1; %nó central
        inci(i,2) = i+1; %primeiro nó
        inci(i,3) = i+2; %segundo nó
    end
    
    figure(6)
    if nnos == 11
         for i = 1:size(inci,1)
             x = [coord(inci(i,1),2),coord(inci(i,2),2),coord(inci(i,3),2)];
             y = [coord(inci(i,1),3),coord(inci(i,2),3),coord(inci(i,3),3)];
             if mod(i,2) == 0
                patch(x,y,[0 0 0])
             else
                 patch(x,y,[1 1 1])
             end
         end
    end

    %% Cálculo da área do quarto de círculo

    x1 = coord(inci(1,1),2); %coordenadas do ponto central
    y1 = coord(inci(1,1),3);

    f = zeros(nnos-1,1);

    for i = 1:nnos-1
        x2 = coord(inci(i,2),2); %coordenadas x dos pontos localizados na circunferência
        x3 = coord(inci(i,3),2);

        y2 = coord(inci(i,2),3); %coordenadas y dos pontos localizados na circunferência
        y3 = coord(inci(i,3),3);

        a = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)); %cálculo da área do elemento

        f(i) = a;

        A = A + a; %somatório da área dos elementos
    end

    erro_Area = norm(((pi)/4)-A);

    %% Calculo do Centróide 

    k = 0;

    j = 0;

    for i = 1:nnos-1
        x2 = coord(inci(i,2),2); %coordenadas x dos pontos localizados na circunferência
        x3 = coord(inci(i,3),2);

        y2 = coord(inci(i,2),3); %coordenadas y dos pontos localizados na circunferência
        y3 = coord(inci(i,3),3);

        k = k + ((x1+x2+x3)/3)*f(i);

        j = j + ((y1+y2+y3)/3)*f(i); 
    end

    X = k/A; %centróide X

    Y = j/A; %centróide Y

    Centroide_exato = (4/(pi*3))+coord(1,2);

    erro_X = norm(Centroide_exato-X); %erro do calculo do valor de X do centroide
    erro_Y = norm(Centroide_exato-Y); %erro do calculo do valor de Y do centroide

    %% Calculo do momento de inércia do quarto de circulo

    Ix = 0; %momento de inércia em relação ao eixo x

    Iy = 0; %momento de inércia em relação ao eixo y

    for i = 1:nnos-1
        x2 = coord(inci(i,2),2); %coordenadas x dos pontos localizados na circunferência
        x3 = coord(inci(i,3),2);

        y2 = coord(inci(i,2),3); %coordenadas y dos pontos localizados na circunferência
        y3 = coord(inci(i,3),3);

        Ix = Ix + 1/12*(y1^2+y2^2+y3^2+y1*y2+y1*y3+y2*y3)*(x2*y3-x3*y2); %momento de inércia em relação ao eixo x
        Iy = Iy + 1/12*(x1^2+x1*x2+x1*x3+x2^2+x2*x3+x3^2)*(x2*y3-x3*y2); %momento de inércia em relação ao eixo y
    end

    Iexato = pi/16; %valor exato do momento de inércia do quarto de círculo para R = 1

    erro_Ix = norm(Iexato-Ix);
    erro_Iy = norm(Iexato-Iy);
    
    if nnos+1 == 10
        e_Area(1,1) = erro_Area; 
        e_Centroide(1,1) = erro_X;
        e_Inercia(1,1) = erro_Ix;
    elseif nnos+1 == 20
        e_Area(2,1) = erro_Area; 
        e_Centroide(2,1) = erro_X;
        e_Inercia(2,1) = erro_Ix;
    elseif nnos+1 == 30
        e_Area(3,1) = erro_Area; 
        e_Centroide(3,1) = erro_X;
        e_Inercia(3,1) = erro_Ix;
    elseif nnos+1 == 40
        e_Area(4,1) = erro_Area; 
        e_Centroide(4,1) = erro_X;
        e_Inercia(4,1) = erro_Ix;
    elseif nnos+1 == 50
        e_Area(5,1) = erro_Area; 
        e_Centroide(5,1) = erro_X;
        e_Inercia(5,1) = erro_Ix;
    elseif nnos+1 == 55
        e_Area(6,1) = erro_Area; 
        e_Centroide(6,1) = erro_X;
        e_Inercia(6,1) = erro_Ix;
    end
    nnos = nnos + 1;
end

figure(1) %janela da figura 1
xlim([0 1])
ylim([0 1])
scatter(coord(:,2),coord(:,3),50,[1 0 0]) %gráfico dos nós
title('Nós do quarto de círculo')
xlabel('coordenadas x dos nós')
ylabel('coordenadas y dos nós')


figure(2) %janela da figura 2
fill(coord(:,2),coord(:,3),[0 0 0]) %gráfico de área definida pelos nós
title('Área do quarto de círculo')
xlabel('coordenadas x dos nós')
ylabel('coordenadas y dos nós')

figure(3) %janela da figura 3
x = [10,20,30,40,50,55];
plot(x,e_Area)
title('Erro da área')
xlabel('Quantidade de nós')
ylabel('\phi')

figure(4) %janela da figura 4
x = [10,20,30,40,50,55];
plot(x,e_Centroide)
title('Erro do centróide')
xlabel('Quantidade de nós')
ylabel('\phi')

figure(5) %janela da figura 5
x = [10,20,30,40,50,55];
plot(x,e_Inercia)
title('Erro da inércia')
xlabel('Quantidade de nós')
ylabel('\phi')