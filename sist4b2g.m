%% Dimensionamento do sistema de pot�ncia
%
% n: n�mero de barras; nl: n�mero de ramos (linhas e trafos); ng: n�mero de 
%                                                                 geradores
%
n = 5; nl = 6; ng = 3;
%
%% Dados dos geradores
%  Compostos por:
%  * Barra �'k' onde unidade geradora est� conectada;
%  * Dados da Fun��o de Produ��o F(P)= c0 + c1 * P + (1/2) * c2 * P^2
%    (assume-se que P est� expresso em MW);      
%  * Limites de gera��o: p_min e p_max, ambos em MW.
%  Cada linha da matriz 'dad_ger' corresponde a um gerador. � formada por:
%  k, c0, c1, c2, p_min, p_max
%
dad_ger=[1 9180 129.6 0.0512 150 600;
         3 6200 157 0.0776 100 400;
         4 1560 159.4 0.01928 50 200;
         ];  
%
%% Dados de carga por barra
%  Matriz n x 2, cuja 1a. coluna corresponde ao n�mero identificador das 
%                barras e 2a. coluna � a respectiva carga, em MW
%
dad_bar =[1 0;
          2 300;
          3 150;
          4 200;
          5 200];
%
%% Par�metros dos ramos (linhas e trafos)
%
% Cada linha da matriz 'dad_lin' corresponde a um ramo. � formada por:
% Barra de origem, barra de chegada, resist�ncia (em p.u.), reat�ncia
% indutiva (em p.u.), limites m�nimo e m�ximo de fluxo (em MW)
%
dad_lin =[1 2 0.02 0.20 -700 700;
          2 3 0.03 0.30 -700 700;
          3 5 0.02 0.20 -700 700;
          3 4 0.04 0.40 -700 700;
          5 4 0.063 0.252 -700 700;
          1 5 0.031 0.126 -170 170];
%
%
