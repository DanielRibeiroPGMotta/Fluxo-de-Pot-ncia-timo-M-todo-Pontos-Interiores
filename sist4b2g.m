%% Dimensionamento do sistema de potência
%
% n: número de barras; nl: número de ramos (linhas e trafos); ng: número de 
%                                                                 geradores
%
n = 5; nl = 6; ng = 3;
%
%% Dados dos geradores
%  Compostos por:
%  * Barra ´'k' onde unidade geradora está conectada;
%  * Dados da Função de Produção F(P)= c0 + c1 * P + (1/2) * c2 * P^2
%    (assume-se que P está expresso em MW);      
%  * Limites de geração: p_min e p_max, ambos em MW.
%  Cada linha da matriz 'dad_ger' corresponde a um gerador. É formada por:
%  k, c0, c1, c2, p_min, p_max
%
dad_ger=[1 9180 129.6 0.0512 150 600;
         3 6200 157 0.0776 100 400;
         4 1560 159.4 0.01928 50 200;
         ];  
%
%% Dados de carga por barra
%  Matriz n x 2, cuja 1a. coluna corresponde ao número identificador das 
%                barras e 2a. coluna é a respectiva carga, em MW
%
dad_bar =[1 0;
          2 300;
          3 150;
          4 200;
          5 200];
%
%% Parâmetros dos ramos (linhas e trafos)
%
% Cada linha da matriz 'dad_lin' corresponde a um ramo. É formada por:
% Barra de origem, barra de chegada, resistência (em p.u.), reatância
% indutiva (em p.u.), limites mínimo e máximo de fluxo (em MW)
%
dad_lin =[1 2 0.02 0.20 -700 700;
          2 3 0.03 0.30 -700 700;
          3 5 0.02 0.20 -700 700;
          3 4 0.04 0.40 -700 700;
          5 4 0.063 0.252 -700 700;
          1 5 0.031 0.126 -170 170];
%
%
