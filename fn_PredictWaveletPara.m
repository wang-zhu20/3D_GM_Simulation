%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get predicted wavelet parameter from emperical equation
%%% Sept 21, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Note: totalEnergy should *0.01
function [outprm]=fn_PredictWaveletPara(M,Rrup,Rhyp,Vs30)  %%%% two columns
% M=7;
% Rrup=100;
% Rhyp=100;
% Vs30=270;
outprm = fn_initDBs(3);
%%
load covrho_3dim.mat

%%%%%%%%%%%%%%%%%%  Baker's Prediction Equation          %%%%%%%%%%%%%%%%
%%%% Table 3.2: Coefficients of the prediction equation, p. 59  %%%%%%%%%
%%%% alpha  beta1  beta2  beta3  beta4  beta5  beta6  h  sigma  tau  %%%%


flag=1;

%%%% WZ's regression %%%%%%%%%%%%%%%%%
while true
    
    epst = mvnrnd(zeros(36,1),cov_total);
    
    coeff = [2.32458719	0	0	0.000506597	-0.000539668	0.30642475	-0.201291382	1	0.2400823	0.2778728
        2.580219233	0	0	0.000514262	-0.00657608	0.193862633	-0.202233797	1	0.3458085	0.3085608
        -0.644272974	0.039631842	0	0	-0.008916714	-0.090199703	0.419026062	10	0.383876	0.2838455
        -0.581144738	0.176824238	0	0	-0.007690044	-0.165130401	0.318061757	10	0.4941612	0.3414261
        -0.594955203	0.036658372	0	0	-0.000829481	-0.043644695	0.063148914	10	0.06321944	0.04280378
        2.049512881	0	0	0.000540427	-0.000302362	0.362411102	-0.218692462	1	0.2575868	0.294046
        2.456448711	0	0	0.000549458	-0.006366325	0.239812844	-0.317116381	1	0.4023611	0.2903117
        -0.886265788	-0.049173756	0	0	-0.008897877	-0.026121137	0.461098691	10	0.4093361	0.2807718
        -0.750156155	0.000106757	0	0	-0.007776414	-0.082442628	0.336617546	10	0.4443033	0.2868338
        -1.054990105	0.033012506	0	0	-0.000271073	-0.054167138	0.129906626	10	0.1473161	0.06819508
        -44.82291085	-5.744236114	45.88888037	0	0	-2.237891286	-0.92030289	10	1.222418	0.6286348
        -36.07664957	-4.241656818	37.73739609	0	0	-2.040610427	-0.736330486	10	0.9269778	0.5765289
        1.2602	0	0	0	0	0	0	0	0.1549	0];

    alpha=coeff(:,1);
    beta1=coeff(:,2);
    beta2=coeff(:,3);
    beta3=coeff(:,4);
    beta4=coeff(:,5);
    beta5=coeff(:,6);
    beta6=coeff(:,7);
    h=coeff(:,8);
    sigma=coeff(:,9);
    tau=coeff(:,10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      for jj=1:24
    % intra_w(jj)=sqrt(covm(jj,jj));
    % inter_w(jj)=sqrt(covi(jj,jj));
    %      end
    %%%%%%%%%%%%%%%%%%%  Compute Median Prediction of Wavelet Parameters %%%
    
    Y=alpha+beta1*M + beta2*log(M)+beta3*exp(M)+beta4*(Rhyp-Rrup)+beta5.*log(sqrt(Rrup.^2+h.^2))+beta6*log(Vs30);
    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12));
        majorEa  = exp(Y(11));
    else
        
        minorEx = exp(Y(1)+epst(1));
        minorSx = exp(Y(2)+epst(2));
        minorEy = exp(Y(3)+epst(3));
        minorSy = exp(Y(4)+epst(4));
        minorRxy = Y(5)+epst(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(6));
        majorSx = exp(Y(7)+epst(7));
        majorEy = exp(Y(8)+epst(8));
        majorSy = exp(Y(9)+epst(9));
        majorRxy =Y(10)+epst(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12)+epst(12));
        majorEa  = exp(Y(11)+epst(11));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=coeff(13,1);
    sn13=coeff(13,9);
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err ~= 0
        continue
    end
    outprm(1) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
        'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
        'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
        'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd);
    %%
    
    Y=alpha+beta1*M + beta2*log(M)+beta3*exp(M)+beta4*(Rhyp-Rrup)+beta5.*log(sqrt(Rrup.^2+h.^2))+beta6*log(Vs30);
    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12));
        majorEa  = exp(Y(11));
    else
        
        minorEx = exp(Y(1)+epst(1+12));
        minorSx = exp(Y(2)+epst(2+12));
        minorEy = exp(Y(3)+epst(3+12));
        minorSy = exp(Y(4)+epst(4+12));
        minorRxy = Y(5)+epst(5+12);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(6+12));
        majorSx = exp(Y(7)+epst(7+12));
        majorEy = exp(Y(8)+epst(8+12));
        majorSy = exp(Y(9)+epst(9+12));
        majorRxy =Y(10)+epst(10+12);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12)+epst(12+12));
        majorEa  = exp(Y(11)+epst(11+12));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=coeff(13,1);
    sn13=coeff(13,9);
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err ~= 0
        continue
    end
    outprm(2) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
        'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
        'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
        'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd);
    
    %%  for vertical
    coeff = [1.822627315	0	0	0.000516032	-0.001061341	0.311184431	-0.134290412	1	0.2445434	0.2654552
        2.330153913	0	0	0.000498358	-0.005818246	0.306695394	-0.205655346	1	0.3344823	0.3145754
        2.033821664	0.118103742	0	0	-0.008915903	-0.277002987	0.081522062	10	0.3705392	0.2889316
        2.189469896	0.236365575	0	0	-0.007116459	-0.326874568	-0.007297691	10	0.461751	0.3139258
        -0.475909042	0.040868584	0	0	-0.00108546	-0.025588532	0.027868227	10	0.04399827	0.03726009
        1.51937839	0	0	0.000543657	6.60078E-05	0.357584725	-0.148099577	1	0.2663513	0.2753488
        2.151838245	0	0	0.000475269	-0.002859721	0.311682418	-0.292280378	1	0.35067	0.2793521
        1.946551276	0.065582506	0	0	-0.010284787	-0.283370829	0.100615346	10	0.4174516	0.3220831
        1.841882796	0.075719005	0	0	-0.008342107	-0.285399167	0.036062152	10	0.4197931	0.2892622
        -0.9671687	0.058205581	0	0	-0.001617681	-0.054923755	0.080195344	10	0.1175948	0.08218528
        -49.59771459	-5.91496767	46.61493695	0	0	-2.322809902	-0.374791801	10	1.026542	0.5095485
        -36.82097327	-3.904360626	35.90187068	0	0	-2.242173903	-0.454544944	10	0.8420813	0.5974369
        1.3361	0	0	0	0	0	0	0	0.1837	0];
    alpha=coeff(:,1);
    beta1=coeff(:,2);
    beta2=coeff(:,3);
    beta3=coeff(:,4);
    beta4=coeff(:,5);
    beta5=coeff(:,6);
    beta6=coeff(:,7);
    h=coeff(:,8);
    sigma=coeff(:,9);
    tau=coeff(:,10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%  Compute Median Prediction of Wavelet Parameters %%%
    
    Y=alpha+beta1*M + beta2*log(M)+beta3*exp(M)+beta4*(Rhyp-Rrup)+beta5.*log(sqrt(Rrup.^2+h.^2))+beta6*log(Vs30);
    if flag==0
        minorEx = exp(Y(1));
        minorSx = exp(Y(2));
        minorEy = exp(Y(3));
        minorSy = exp(Y(4));
        minorRxy = Y(5);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6));
        majorSx = exp(Y(7));
        majorEy = exp(Y(8));
        majorSy = exp(Y(9));
        majorRxy =Y(10);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12));
        majorEa  = exp(Y(11));
    else
        minorEx = exp(Y(1)+epst(13+12));
        minorSx = exp(Y(2)+epst(14+12));
        minorEy = exp(Y(3)+epst(15+12));
        minorSy = exp(Y(4)+epst(16+12));
        minorRxy = Y(5)+epst(17+12);      minorRxy =2*normcdf(minorRxy,0,1)-1;
        majorEx = exp(Y(6)+epst(18+12));
        majorSx = exp(Y(7)+epst(19+12));
        majorEy = exp(Y(8)+epst(20+12));
        majorSy = exp(Y(9)+epst(21+12));
        majorRxy =Y(10)+epst(22+12);        majorRxy=2*normcdf(majorRxy,0,1)-1;
        totalEnergy = exp(Y(12)+epst(24+12));
        majorEa  = exp(Y(11)+epst(23+12));
    end
    %%%% Version 1: random is longmornal
    %minorRnd = exp(Y(13));
    mn13=coeff(13,1);
    sn13=coeff(13,9);
    m13=log(mn13^2/sqrt(sn13^2+mn13^2));
    s13=sqrt(log(sn13^2/mn13^2+1));
    minorRnd=lognrnd(m13,s13);
    
    %% to lxly
    majorVx = majorSx^2;
    majorVy = majorSy^2;
    majorExy = majorRxy*majorSx*majorSy + majorEx*majorEy;
    
    majorElx = log(majorEx^2 / sqrt(majorVx+majorEx^2));
    majorVlx = log(majorVx/majorEx^2 + 1);
    majorSlx = sqrt(majorVlx);
    majorEly = log(majorEy^2 / sqrt(majorVy+majorEy^2));
    majorVly = log(majorVy/majorEy^2 + 1);
    majorSly = sqrt(majorVly);
    
    majorCovlxly = log(majorExy/majorEx/majorEy);
    majorRlxly = majorCovlxly/majorSlx/majorSly;
    
    minorVx = minorSx^2;
    minorVy = minorSy^2;
    minorExy = minorRxy*minorSx*minorSy + minorEx*minorEy;
    
    minorElx = log(minorEx^2 / sqrt(minorVx+minorEx^2));
    minorSlx = sqrt(log(minorVx/minorEx^2 + 1));
    minorVlx = minorSlx^2;
    minorEly = log(minorEy^2 / sqrt(minorVy+minorEy^2));
    minorSly = sqrt(log(minorVy/minorEy^2 + 1));
    minorVly = minorSly^2;
    
    minorCovlxly = log(minorExy/minorEx/minorEy);
    minorRlxly = minorCovlxly/minorSlx/minorSly;
    
    majorLLCov=[majorVlx majorCovlxly; majorCovlxly majorVly];
    [T,err] = cholcov(majorLLCov);
    if err == 0
        break
    end
end
outprm(3) = struct( 'M',M,'hdist',Rrup,'vs30',Vs30,...
    'minorElx',minorElx, 'minorSlx',minorSlx, 'minorVlx',minorVlx, 'minorEly',minorEly, 'minorSly',minorSly, 'minorVly',minorVly,'minorRlxly',minorRlxly, ...
    'majorElx',majorElx, 'majorSlx',majorSlx, 'majorVlx',majorVlx, 'majorEly',majorEly, 'majorSly',majorSly, 'majorVly',majorVly,'majorRlxly',majorRlxly, ...
    'totalEnergy',totalEnergy, 'majorEa',majorEa, 'minorRnd',minorRnd);

%% save('MedianWaveletParam.mat','MedianWaveletParam');
