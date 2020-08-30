%% Program: Análise térmica 
%
%   @DESCRIÇÃO Atividade 7: Realizar a análise térmica por elemento
%   triangular
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar memória e fechar janelas
clear 
close all
clc

%% Matriz de Coordenadas do sistema 
%   coord=[nº | X | Y] => Matriz de coordenadas dos nós

coord = zeros(36,3);
k = 1;

for i = 0:1:5 %coordenada y 
    for j = 0:1:10 %coordenada x
        coord(k,:) = [k i j];
        k = k + 1; %quantidade de nós
    end
end

nnos = size(coord,1); %número de nós da matriz de coordenadas

l = abs(coord(1,3)-coord(2,3)); %distância nodal

coord2 = zeros(4,3);
k = 1;

for i = 0:1:2 %coordenada y 
    for j = 3:1:7 %coordenada x
        coord2(k,:) = [k i j];
        k = k + 1; %quantidade de nós
    end
end

nnos1 = size(coord,1)-size(coord2,1);
nnos2 = size(coord2,1);

for i = 1 : nnos1
    for j = 1 : nnos2
        if coord(i,2:3) == coord2(j,2:3)
            coord(i,:) = [];
        end
    end
end

centro = [0 5];

pontos = 30;

theta = (pi)/(pontos-1);

nnos3 = size(coord,1);

for i = 1 : pontos
coord(nnos3+i,1) = nnos+i; %número do nó
coord(nnos3+i,3) = 2*cos(theta*(i-1)) + centro(2); %coordenada y
coord(nnos3+i,2) = 2*sin(theta*(i-1)) + centro(1); %coordenada x
end

for i = 1 : size(coord,1)
    coord(i,1) = i;
end

% scatter(coord(:,3),coord(:,2),7,'b')

 condcont = (coord(:,3) == 0); 
 bc = coord(condcont>0,:); %condições de contorno

%% Matriz de incidência
% inci=[nº | tab.mat | tab.geo | nó1 | nó2 ] => Matriz de incidência do elemento

inci = delaunay(coord(:,3),coord(:,2));

apagar=zeros(size(inci,1),1); 

for ii = 1:size(inci,1)
      no1 = inci(ii,1); %nós do elemento triangular
      no2 = inci(ii,2);
      no3 = inci(ii,3);
      
      x1 = coord(no1,3); %coordenada x do nó 1
      x2 = coord(no2,3); %coordenada x do nó 2
      x3 = coord(no3,3); %coordenada x do nó 3
    
      y1 = coord(no1,2); %coordenada y do nó 1
      y2 = coord(no2,2); %coordenada y do nó 2
      y3 = coord(no3,2); %coordenada y do nó 3
      
      x = (x1+x2+x3)/3;
      y = (y1+y2+y3)/3;
      if x^2+y^2 <= 49 && x^2+y^2 >= 9
          if y < 2 && x > 3 &&  x < 7
              apagar(ii,:)=500;
          end
      end
end

inci(apagar>499,:)=[];

for i = 1:size(inci,1) %patch dos elementos
    x = [coord(inci(i,1),3),coord(inci(i,2),3),coord(inci(i,3),3)];
    y = [coord(inci(i,1),2),coord(inci(i,2),2),coord(inci(i,3),2)];

    patch(x,y,[1 1 1])
end

chr = int2str(coord(:,1));
text(coord(:,3),coord(:,2),chr)

%%  condições de contorno
%
%   bc = [nó | grau de liberdade | valor]
%
%   Grau de liberdade 1 --> x
%   Grau de liberdade 2 --> y
%   Grau de liberdade 3 --> z

for i = 1 : size(bc)
    bc(i,2) = 1;
    bc(i,3) = 0;
end

%% fluxo do contorno

condcont = (coord(:,3) == 10);
ccont = coord(condcont>0,:);

%% Matriz de identificação

ngdl = 1;
nnos = size(coord,1);

id = ones(ngdl,nnos);
nbc = size(bc,1); %numero de condições de contorno

for i=1:nbc
    id(bc(i,2),bc(i,1)) = bc(i,3);
end

neq = 0;

for i= 1:nnos
    for j = 1:ngdl
        if id(j,i)== 1
            neq = neq +1;
            id(j,i) = neq;
        end
    end
end

%% Calculo das matrizes e vetores

      q = 50; %carga nodal
      q2 = 10;%carga nodal
      fg = zeros(neq,1); % pré locação do vetor de cargas nodais
      fgc = zeros(neq,1); % pré locação do vetor de cargas nodais do contorno
      nel = size(inci,1);
      kg = zeros(neq,neq); % pré locação da matriz de rigidez global
      
     for ii = 1 : nel
          
      no1 = inci(ii,1); %nós do elemento triangular
      no2 = inci(ii,2);
      no3 = inci(ii,3);
      
      x1 = coord(no1,3); %coordenada x do nó 1
      x2 = coord(no2,3); %coordenada x do nó 2
      x3 = coord(no3,3); %coordenada x do nó 3
    
      y1 = coord(no1,2); %coordenada y do nó 1
      y2 = coord(no2,2); %coordenada y do nó 2
      y3 = coord(no3,2); %coordenada y do nó 3
      
      x = [x1,x2,x3];
      y = [y1,y2,y3];
      a = polyarea(x,y); %área dos elementos
      
      P1 = [x1,y1]; %localização dos nós
      P2 = [x2,y2];
      P3 = [x3,y3];
      
      beta1 = atan2(2*a,dot(P2-P1,P3-P1)); %ângulo do nó 1
      beta2 = atan2(2*a,dot(P1-P2,P3-P2)); %ângulo do nó 2
      beta3 = pi - abs(beta1) - abs(beta2);%ângulo do nó 3

%% matriz de rigidez do elemento 
    
    i = cot(abs(beta1)); %ângulos do elemento
    j = cot(abs(beta2));
    k = cot(abs(beta3));
    
    ke = 0.5*[j+k, -k, -j; %matriz de rigidez local
               -k, i+k, i;
                -j, -i, i+j];    
            
%% matriz de rigidez global  
            
    loc=[id(1,no1) id(1,no2) id(1,no3)]; %matriz de localização da matriz de rigidez local na global
      
    ngdl = 3;
    
    for j = 1:ngdl
        if loc(j) ~= 0
            for k = 1:ngdl
                if loc(k) ~=0
                    kg(loc(j),loc(k))=kg(loc(j),loc(k))+ke(j,k); %matriz de rigidez global
                end
            end
        end
    end    
    
 %% vetor de cargas nodais global   
    
     fe = q*a/3*[1;1;1]; %vetor de cargas nodais local
    
    for jj = 1:ngdl
        if loc(jj) ~= 0
            fg(loc(jj),1)=fg(loc(jj),1)+fe(jj,1); %vetor de cargas nodais global
        end
    end
    
    %% vetor de cargas nodais global contorno
    
    fc = q2*l/2*[1;1];
    
    for kk = 1:size(ccont,1)
        for jj = 1 : size(ccont,1)
            if no1 == ccont(kk,1) && no2 == ccont(jj,1)
                loc = [id(1,no1) id(1,no2)];
                for i = 1:2
                   if loc(i) ~= 0
                     fgc(loc(i),1)=fgc(loc(i),1)+fc(i,1); %vetor de cargas nodais do contorno
                   end
                end 
            elseif no1 == ccont(kk,1) && no3 == ccont(jj,1)
                loc = [id(1,no1) id(1,no3)];
                for i = 1:2
                   if loc(i) ~= 0
                     fgc(loc(i),1)=fgc(loc(i),1)+fc(i,1); %vetor de cargas nodais do contorno
                   end
                end
            elseif no2 == ccont(kk,1) && no3 == ccont(jj,1)
                loc = [id(1,no2) id(1,no3)];
                for i = 1:2
                   if loc(i) ~= 0
                     fgc(loc(i),1)=fgc(loc(i),1)+fc(i,1); %vetor de cargas nodais do contorno
                   end
                end
            end
        end
    end
    
     end
      
     F = fgc + fg; %carga nodal total
     
     T = zeros(size(coord,1),1);
    for i = 1:size(coord,1)
        if i > nnos3
            T(i,1) = 100;
        end
    end
    
    apagar = (id==0);
    T(apagar>0,:)=[]; %eliminando linhas do vetor temperatura
    
    F=F-kg*T; %carga nodal 
    
    remove = (T>99); %remoção das linhas de cargas não nulas
    
    F(remove>0,:) = [];
    kg(remove>0,:) = [];
    kg(:,remove>0) = [];

    T = kg\F;
    
    recuperar = (id==0);
    T0 = T;
    T0(recuperar>0) = 0;
    
       
     for i = 1:size(coord,1)
        if i > nnos3
            T0(i,1) = 100;
        end
     end
    
     figure(2)
     scatter(coord(:,3),coord(:,2),7,'b')
     chr = int2str(T0(:,1));
     text(coord(:,3),coord(:,2),chr)
     
     %         if no1 == ccont(kk,1)
%             for jj = 1:3
%                 if no2 == ccont(jj,1)
%                     loc = [id(1,no1) id(1,no2)];
%                     for i = 1:2
%                        if loc(i) ~= 0
%                          fgc(loc(i),1)=fgc(loc(i),1)+fc(i,1); %vetor de cargas nodais do contorno
%                        end
%                     end                    
%                 elseif no3 == ccont(jj,1)
%                     loc = [id(1,no1) id(1,no3)];
%                     for i = 1:2
%                        if loc(i) ~= 0
%                          fgc(loc(i),1)=fgc(loc(i),1)+fc(i,1); %vetor de cargas nodais do contorno
%                        end
%                     end
%                 end
%             end
%         elseif no2 == ccont(kk,1)
%             for jj = 1:3
%                 if no3 == ccont(jj,1)
%                     loc = [id(1,no2) id(1,no3)];
%                     for i = 1:2
%                        if loc(i) ~= 0
%                          fgc(loc(i),1)=fgc(loc(i),1)+fc(i,1); %vetor de cargas nodais do contorno
%                        end
%                     end
%                 end
%             end
%         end
