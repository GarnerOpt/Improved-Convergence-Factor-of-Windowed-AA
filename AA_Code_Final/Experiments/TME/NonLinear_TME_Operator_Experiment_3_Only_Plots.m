%==========================================================================
% TME Experiment 3 (Plots Only)
%==========================================================================

load('TME_DM1_p100_n105')


%==========================================================================
% Make Box-Plot of Numerial Experiments
%==========================================================================
indx=find(Methods==1);
Times_AA = Times;
Times_AA(indx)=[];

Iterations_AA = Iterations;
Iterations_AA(indx)=[];

Methods_AA = Methods;
Methods_AA(indx)=[];

Alg_AA = categorical(Methods_AA,[2,3,4],{'AA(1)','AA(2)','AA(3)'});
Alg = categorical(Methods,[1,2,3,4],{'FP','AA(1)','AA(2)','AA(3)'});

%--------------------------------------------------------------------------
% Boxplot for Compute Time
%--------------------------------------------------------------------------
figure;
b1 = axes();
boxplot(b1,Times,Alg)
title('Compute Time Comparison','FontSize',16,'FontWeight','bold')
xlabel('Method','FontSize',16,'FontWeight','bold')
ylabel('Time (s)','FontSize',16,'FontWeight','bold')
ax3 = gca;
ax3.FontWeight = 'bold';


%--------------------------------------------------------------------------
% Boxplot for Total Iterations
%--------------------------------------------------------------------------
figure;
b1 = axes();
boxplot(b1,Iterations,Alg)
title('Total Iteration Comparison','FontSize',16,'FontWeight','bold')
xlabel('Method','FontSize',16,'FontWeight','bold')
ylabel('Iterations','FontSize',16,'FontWeight','bold')
ax1 = gca;
ax1.FontWeight = 'bold';

