%% Program: �rea, Centr�ide e Momento de in�rcia
%
%   @DESCRI��O Atividade 5: Definir a �rea da geometria, seu centr�ide e
%   momento de in�rcia
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar mem�ria e fechar janelas
clear all
close all
clc

%% Matriz de Coordenadas do sistema 
%   coord=[n� | X | Y] => Matriz de coordenadas dos n�s

nnos = 5; %n�mero de n�s menos 1
nel = nnos - 1;

A = 0; %somat�rio das �reas dos elementos
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
    coord(i+1,1) = i+1; %n�mero do n�
    coord(i+1,2) = cos(theta*(i-1))+coord(1,2); %coordenada x
    coord(i+1,3) = sin(theta*(i-1))+coord(1,3); %coordenada y
    end



    %-------Matriz de incid�ncia ---------------------------------------------%
    % inci=[n� | tab.mat | tab.geo | n�1 | n�2 ] => Matriz de
    % incid�ncia do elemento

    inci = zeros(nnos-1,3); 

    for i = 1:nnos-1 
        inci(i,1) = 1; %n� central
        inci(i,2) = i+1; %primeiro n�
        inci(i,3) = i+2; %segundo n�
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

    %% C�lculo da �rea do quarto de c�rculo

    x1 = coord(inci(1,1),2); %coordenadas do ponto central
    y1 = coord(inci(1,1),3);

    f = zeros(nnos-1,1);

    for i = 1:nnos-1
        x2 = coord(inci(i,2),2); %coordenadas x dos pontos localizados na circunfer�ncia
        x3 = coord(inci(i,3),2);

        y2 = coord(inci(i,2),3); %coordenadas y dos pontos localizados na circunfer�ncia
        y3 = coord(inci(i,3),3);

        a = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)); %c�lculo da �rea do elemento

        f(i) = a;

        A = A + a; %somat�rio da �rea dos elementos
    end

    erro_Area = norm(((pi)/4)-A);

    %% Calculo do Centr�ide 

    k = 0;

    j = 0;

    for i = 1:nnos-1
        x2 = coord(inci(i,2),2); %coordenadas x dos pontos localizados na circunfer�ncia
        x3 = coord(inci(i,3),2);

        y2 = coord(inci(i,2),3); %coordenadas y dos pontos localizados na circunfer�ncia
        y3 = coord(inci(i,3),3);

        k = k + ((x1+x2+x3)/3)*f(i);

        j = j + ((y1+y2+y3)/3)*f(i); 
    end

    X = k/A; %centr�ide X

    Y = j/A; %centr�ide Y

    Centroide_exato = (4/(pi*3))+coord(1,2);

    erro_X = norm(Centroide_exato-X); %erro do calculo do valor de X do centroide
    erro_Y = norm(Centroide_exato-Y); %erro do calculo do valor de Y do centroide

    %% Calculo do momento de in�rcia do quarto de circulo

    Ix = 0; %momento de in�rcia em rela��o ao eixo x

    Iy = 0; %momento de in�rcia em rela��o ao eixo y

    for i = 1:nnos-1
        x2 = coord(inci(i,2),2); %coordenadas x dos pontos localizados na circunfer�ncia
        x3 = coord(inci(i,3),2);

        y2 = coord(inci(i,2),3); %coordenadas y dos pontos localizados na circunfer�ncia
        y3 = coord(inci(i,3),3);

        Ix = Ix + 1/12*(y1^2+y2^2+y3^2+y1*y2+y1*y3+y2*y3)*(x2*y3-x3*y2); %momento de in�rcia em rela��o ao eixo x
        Iy = Iy + 1/12*(x1^2+x1*x2+x1*x3+x2^2+x2*x3+x3^2)*(x2*y3-x3*y2); %momento de in�rcia em rela��o ao eixo y
    end

    Iexato = pi/16; %valor exato do momento de in�rcia do quarto de c�rculo para R = 1

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
scatter(coord(:,2),coord(:,3),50,[1 0 0]) %gr�fico dos n�s
title('N�s do quarto de c�rculo')
xlabel('coordenadas x dos n�s')
ylabel('coordenadas y dos n�s')


figure(2) %janela da figura 2
fill(coord(:,2),coord(:,3),[0 0 0]) %gr�fico de �rea definida pelos n�s
title('�rea do quarto de c�rculo')
xlabel('coordenadas x dos n�s')
ylabel('coordenadas y dos n�s')

figure(3) %janela da figura 3
x = [10,20,30,40,50,55];
plot(x,e_Area)
title('Erro da �rea')
xlabel('Quantidade de n�s')
ylabel('\phi')

figure(4) %janela da figura 4
x = [10,20,30,40,50,55];
plot(x,e_Centroide)
title('Erro do centr�ide')
xlabel('Quantidade de n�s')
ylabel('\phi')

figure(5) %janela da figura 5
x = [10,20,30,40,50,55];
plot(x,e_Inercia)
title('Erro da in�rcia')
xlabel('Quantidade de n�s')
ylabel('\phi')