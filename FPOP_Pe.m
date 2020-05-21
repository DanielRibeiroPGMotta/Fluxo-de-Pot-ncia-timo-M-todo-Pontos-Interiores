% FPOP-PE: solução via método primal-dual de pontos interiores
%          do problema de minimização dos custos de produção
%          para fluxo de potência ativa, com representação da
%          rede e consideração das perdas de transmissão
%          (despacho reativo não é considerado)
% Arquivo de dados: dadN*.m, onde N é o número de barras
% Arquivos script necessários: incidmat.m, slacks.m, admitt.m
% Arquivos Function necessários: flows.m, gtheta.m, gp.m tfun.m, 
%                                ufun.m, hfun.m
%
% entrada de dados (originalmente compatíveis com potências em MW)
%
sist4b2g
%dad5bmod
%Dad5b_g
%resistsist5b
%sist5bww
%Dad5b_3g_wwB_MW
%exerc_comp4
%
% Recuperação dos dados de barra
%
ip = dad_ger(:,1);
c0= dad_ger(:,2);
c = dad_ger(:,3);
q = dad_ger(:,4);
pubar = dad_ger(:,5);
pbar = dad_ger(:,6);
%
% Recuperação dos dados de carga
%
d=dad_bar(:,2);
%
% Recuperação dos dados dos ramos
%
na = dad_lin(:,1);
nb = dad_lin(:,2);
RL = dad_lin(:,3);
XL = dad_lin(:,4);
fubar = dad_lin(:,5);
fbar = dad_lin(:,6);
%
%    Conversão de dados para potência em p.u.
%
Sb = 100;
c_0 = c; q_0 = q;
c = c * Sb; q = q * Sb^2;
pubar_0 = pubar; pbar_0 = pbar;
pubar = pubar / Sb; pbar = pbar / Sb;
%
d_0 = d;
d = d / Sb;
%
fubar_0 = fubar; fbar_0 = fbar;
fubar = fubar / Sb; fbar = fbar / Sb;
% p = p / Sb; f = f / Sb;
%
% Inicialização de variáveis
%
theta = zeros(n,1);
p = mean([pubar'; pbar'])';
f = mean([fubar'; fbar'])';
%
num_bar=1:n;
dad_ger=[ip,c0,c_0,q_0,pubar_0,pbar_0];
Q=diag(q);
b = ones(nl,1)./XL;
f=[f;-f];
fbar=[fbar;fbar];
fubar=[fubar;fubar];
Na = [na;nb]; Nb = [nb;na];
%
Incidmat
admitt
% Compute initial values for slacks
Slacks
%
%stop
% Initialize Interior point parameters
mu = 5; sigma = 0.9995; beta = 10; lim_iter = 50; norma = 9999; tol= 1.e-6;
% Initialize Lagrange multipliers for equality constraints
lambda = ones(n,1);
% Compute Lagrange multipliers for inequality constraints
piubar = - mu\subar; 
pibar = mu\sbar;
wubar = - mu\rubar; 
wbar = mu\rbar;
t1 = 1; u1 = 1;
% Initialize counter
counter = 0;
%
% Start main loop
%
while (norma > tol & counter < lim_iter )
%
      counter = counter + 1;
      [f,F2] = Flows (n,nl,na,nb,b,theta,yl_mag,ang_y);
%
      Gp = GP(Q,subar,sbar,piubar,pibar);
      Gtheta = GTHETA(na,nb,b,theta,yl_mag,ang_y,A2,lambda,wubar,wbar,F2,rubar,rbar);
%     
      vubar = mu * ones(ng,1) + subar .* piubar;
      vbar = mu * ones(ng,1) - sbar .* pibar;
      yubar = mu * ones(2*nl,1) + rubar .* wubar;
      ybar = mu * ones(2*nl,1) - rbar .* wbar;
      zubar = f - rubar - fubar;
      zbar = f + rbar -fbar;
%
      [t,t1] = Tfun (c,Q,p,Ag,lambda,piubar,pibar,subar,sbar,vubar,vbar);
      [u,u1] = Ufun (F2,A2,lambda,wubar,wbar,rubar,rbar,yubar,ybar,zubar,zbar);
      h = Hfun (A2,Ag,p,f,d);
%
      AF = A2 * F2;
      M = [Gp zeros(ng,n-1) -Ag';zeros(n-1,ng) Gtheta AF'; 	
           -Ag AF zeros(n,n)];
%
      rhs = [t; u; -h];
%
      x = M\rhs;
      dp = x(1:ng);
      dtheta = x(ng+1:ng+n-1);
      dlambda = x(ng+n:ng+2*n-1);
%
      dsubar = dp;
      dsbar = - dp;
      drubar = F2 * dtheta + zubar;
      drbar = - F2 * dtheta - zbar;
      dpiubar = - subar .\ (piubar .* dp + vubar);
      dpibar = sbar .\ (pibar .* dp + vbar);
      dwubar = - rubar .\ (wubar .* (F2 * dtheta + zubar) + yubar);
      dwbar = rbar .\ (wbar .* (F2 * dtheta + zbar) + ybar);
%
% Computing alphap and alphad
%
      k = find(dsubar<0);
      if isempty(k) == 1
         minsubar = 1;
      else
         ratio = abs(dsubar(k)) .\ subar(k);
         minsubar = min(ratio);
      end
%
      k = find(dsbar<0);
      if isempty(k) == 1
         minsbar = 1;
      else
         ratio = abs(dsbar(k)) .\ sbar(k);
         minsbar = min(ratio);
      end
%      
      k = find(drubar<0);
      if isempty(k) == 1
         minrubar = 1;
      else
         ratio = abs(drubar(k)) .\ rubar(k);
         minrubar = min(ratio);
      end
%	
      k = find(drbar<0);
      if isempty(k) == 1
         minrbar = 1;
      else
         ratio = abs(drbar(k)) .\ rbar(k);
         minrbar = min(ratio);
      end
%
      argp = [minsubar,minsbar,minrubar,minrbar,1];     
      alphap = min (argp);
%
      k = find(dpiubar>0);
      if isempty(k) == 1
         minpiubar = 1;
      else
         ratio = abs(dpiubar(k)) .\ (-piubar(k));
         minpiubar = min(ratio);
      end
%
      k = find(dpibar<0);
      if isempty(k) == 1
         minpibar = 1;
      else
         ratio = abs(dpibar(k)) .\ pibar(k);
         minpibar = min(ratio);
      end
%
      k = find(dwubar>0);
      if isempty(k) == 1
         minwubar = 1;
      else
         ratio = abs(dwubar(k)) .\ (-wubar(k));
         minwubar = min(ratio);
      end
% 
      k = find(dwbar<0);
      if isempty(k) == 1
         minwbar = 1;
      else
         ratio = abs(dwbar(k)) .\ wbar(k);
         minwbar = min(ratio);
      end
%
     argd = [minpiubar,minpibar,minwubar,minwbar,1];     
     alphad = min (argd);
%
% Updating current solution
%
      sap = sigma * alphap;
      sad = sigma * alphad;
%
      p = p + sap * dp;
      theta(2:n) = theta(2:n) + sap * dtheta;
      lambda = lambda + sad * dlambda;
      piubar = piubar + sad * dpiubar;
      pibar = pibar + sad * dpibar;
      wubar = wubar + sad * dwubar;
      wbar = wbar + sad * dwbar;
      subar = subar + sap * dsubar;
      sbar = sbar + sap * dsbar;
      rubar = rubar + sap * drubar;
      rbar = rbar + sap * drbar;
%
% Updating mu
%
	   [f,F] = Flows (n,nl,na,nb,b,theta,yl_mag,ang_y); 
%
      mu_num1 = [subar;sbar;rubar;rbar];
      mu_num2 = [-piubar;pibar;-wubar;wbar];
      mu_num = mu_num1' * mu_num2;
      mu_den = 2 * beta * (ng + n - 1);
      mu = mu_num / mu_den;
      Kktipnt;
      kkT=[kkt1;mu_num];
      norma = norm(kkT,inf);
      followup=[counter;norma];
      fprintf('Criterio de Convergencia na iter %3.0f = %8.3d\n\n',followup);

%
end
cost = c0'*ones(ng,1) + c'* p +0.5 * p'* Q * p;
perdas = (ones(1,ng) * p - ones(1,n) * d) * Sb;
%
%
%
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%
%                           Saida do problema de FPO 
%
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% 
fprintf('==================================================================\n');
fprintf('         Valor final do parametro de barreira = %8.3d\n\n',mu);
if (norma < 1e-6)
  fprintf('         Convergencia obtida em %3.0f iteracoes. \n\n',counter);
else
  fprintf('         Convergencia nao obtida em %3.0f iteracoes\n\n',lim_iter);
end 
fprintf('=================================================================\n');
fprintf('\n\n');

fprintf('=================================================================\n');
%fprintf('-----------------------------------------------------------------\n');
fprintf('                  Componentes da funcao custo               \n');
fprintf('\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('barra   valor constante      coef. linear        coef. quadratico \n'); 
fprintf('-----------------------------------------------------------------\n');
gerout1= (dad_ger(:,1:4))';
gerout = gerout1(:);
fprintf('%3.0f %13.0f %22.4f %20.4f\n',gerout);
disp(' ')
fprintf('=================================================================\n');
%
fprintf('\n\n');
%
fprintf('=======================================================================\n');
fprintf('                        Saída  dos geradores                           \n');
fprintf('=======================================================================\n');
fprintf('                    Geração de Potência Ativa  \n');
fprintf('-----------------------------------------------------------------------\n');
fprintf('barra  Pot. gerada    mult. lagrange    var. folga        limites     \n');
fprintf('          (MW)        sup.    inf.      sup.   inf.      sup.   inf.    \n');
fprintf('-----------------------------------------------------------------------\n');
pbarra=zeros(n,1);
pbarra(ip) = p;
pls=zeros(n,1);
pls(ip)=pbar;
pli=zeros(n,1);
pli(ip)=pubar;
gerativ1=[ ip p*Sb pibar*(1/Sb) piubar*(1/Sb) sbar*Sb subar*Sb pbar_0 pubar_0];
gerativ1=gerativ1';
gerativ=gerativ1(:);
fprintf('%3.0f %12.2f %10.3f %6.3f %9.2f %7.2f %8.2f %6.2f    \n',gerativ);
fprintf('=======================================================================\n');
fprintf('\n\n');
fprintf('=================================================================\n');
disp(' ')
fprintf('            Valor da funcao custo  = %10.2f $/h\n\n',cost);
fprintf('=================================================================\n');
fprintf('\n\n');
fprintf('==================================================================================\n');
fprintf('                             Variaveis das Barras\n');
fprintf('---------------------------------------------------------------------------------\n');
fprintf('barra      tensão          lambda       carga       var. folga (V)    limites (V)    \n');
fprintf('        mag   angulo       ($/MWh)      (MW)        v_sup.  v_inf.   v_sup. v_inf.   \n');
fprintf('---------------------------------------------------------------------------------\n');
j = (1:n)';
ve=ones(n,1);vebar=1.05*ones(n,1);veubar=0.95*ones(n,1);svbar=0.05*ones(n,1);svubar=svbar;
radeg = 180/pi;
barr1=[j ve theta*radeg lambda*(1/Sb) d*Sb svbar svubar vebar veubar];
barr1=barr1';
barr=barr1(:);
fprintf('%3.0f %8.4f %8.2f %11.4f    %8.3f    %9.3f %6.3f %8.2f  %4.2f    \n',barr);
fprintf('==================================================================================\n');
fprintf('\n\n');
[f,F] = Flows (n,nl,na,nb,b,theta,yl_mag,ang_y);
fprintf('==================================================================================\n');
fprintf('                             Variaveis das Linhas\n');
fprintf('----------------------------------------------------------------------------------\n');
fprintf('               Fluxo de pot.       Quantidades referentes ao fluxo ativo          \n');
fprintf('  barras          ativa          multip. lagrange  var. folga       limites   \n');
fprintf('orig. dest.        (MW)           sup.     inf.    sup.    inf.    sup.    inf. \n ');
fprintf('----------------------------------------------------------------------------------\n');
var_lin1=[Na Nb f*Sb wbar*(1/Sb) wubar*(1/Sb) rbar*Sb rubar*Sb fbar*Sb fubar*Sb]; 
var_lin1=var_lin1';
var_lin=var_lin1(:);
fprintf('%4.0d %4.0d     %9.2f      %9.3f %8.3f %8.2f %7.2f  %6.2f %6.2f    \n',var_lin);
fprintf('==================================================================================\n');
disp(' ')
fprintf('                      Perdas de Transmissão = %10.3f MW\n\n',perdas);
fprintf('==================================================================================\n');

% Rotina para Calculo do Encargo e RE
k=1.1;
t=1.2;
Pw=50;

enc21 = k*(lambda(2) - lambda(1))*Pw;
RE21 = t*(lambda(2))*Pw - 0.5*enc21

enc25 = k*(lambda(5) - lambda(1))*Pw;
RE25 = t*(lambda(5))*Pw - 0.5*enc25

enc12 = k*(lambda(2) - lambda(1))*Pw;
RE12 = t*(lambda(2))*Pw - 0.5*enc12

enc13 = k*(lambda(3) - lambda(1))*Pw;
RE13 = t*(lambda(3))*Pw - 0.5*enc13

REbarra1 = RE21 + RE25
REbarra2 = RE12 + RE13
