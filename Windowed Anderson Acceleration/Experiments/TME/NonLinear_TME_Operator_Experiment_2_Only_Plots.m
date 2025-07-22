%==========================================================================
% AISTATS TME Experiment 2 (Plots Only)
%==========================================================================
% GOAL: This script compares AA(m) applied to the alternative TME
% fixed-point iteration with the standard fixed-point iteration for TME for
% Data Model 1. 
%==========================================================================

load('TME_DM2_Exp1')

%======================================
% R-Linear Rate Estimate Plot
%======================================
% Plot Results: R-Linear Convergence Rate
figure(3);
for i=1:100
    plot(R_rate_FP{i},'b-','LineWidth',1.5);
    hold on;
    plot(R_rate_AA{1,i},'r-','LineWidth',1.5);
    hold on;
    plot(R_rate_AA{2,i},'k-','LineWidth',1.5);
    hold on;
    plot(R_rate_AA{3,i},'g-','LineWidth',1.5);
    hold on;
    plot([1,length(Err_FP{1})],[r_est,r_est],'k--','LineWidth',1.5);
    hold on;
end
legend({'FP','AA(1)','AA(2)','AA(3)','Thm 2 Est.'},'FontSize',14,'FontWeight','bold','Location','best');
set(gca, 'YDir', 'reverse')
xlabel('Iterations (k)','FontSize',16,'FontWeight','bold');
ylabel('r_{est}(k)','FontSize',16,'FontWeight','bold')
axis([0,maxAA_iters+10,0.01, 1.5])
title('Comparing R-Linear Convergence Rates','FontSize',16,'FontWeight','bold');


%======================================
% R-Linear Rate Ratio Estimate
%======================================
figure(4);
for i=1:100
    plot([1,length(Err_FP{1})],[RHS,RHS],'k--','LineWidth',2);
    hold on;
    plot(LHS{1,i},'r-','LineWidth',1.5);
    hold on;
    plot(LHS{2,i},'b-','LineWidth',1.5);
    hold on;
    plot(LHS{3,i},'g-','LineWidth',1.5);
    hold on;
end
legend({'RHS','LHS-AA(1)','LHS-AA(2)','LHS-AA(3)'},'FontSize',14,'FontWeight','bold','Location','best');
xlabel('Iterations (k)','FontSize',16,'FontWeight','bold');
ylabel('Values in Eq. (28)','FontSize',16,'FontWeight','bold')
axis([0,maxAA_iters+10,0, 1.1*RHS])
title('Tightness of Ratio','FontSize',16,'FontWeight','bold');


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
figure(5);
b1 = axes();
boxplot(b1,Times,Alg)
title('Compute Time Comparison','FontSize',16,'FontWeight','bold')
xlabel('Method','FontSize',16,'FontWeight','bold')
ylabel('Time (s)','FontSize',16,'FontWeight','bold')
ax3 = gca;
ax3.FontWeight = 'bold';
b2 = axes();
b2.Position = [0.4 0.6 0.5 0.3];
boxplot(b2,Times_AA,Alg_AA)
ax4 = gca;
ax4.FontWeight = 'bold';

%--------------------------------------------------------------------------
% Boxplot for Total Iterations
%--------------------------------------------------------------------------
figure(6);
b1 = axes();
boxplot(b1,Iterations,Alg)
title('Total Iteration Comparison','FontSize',16,'FontWeight','bold')
xlabel('Method','FontSize',16,'FontWeight','bold')
ylabel('Iterations','FontSize',16,'FontWeight','bold')
ax1 = gca;
ax1.FontWeight = 'bold';
b2 = axes();
b2.Position = [0.4 0.6 0.5 0.3];
boxplot(b2,Iterations_AA,Alg_AA)
ax2 = gca;
ax2.FontWeight = 'bold';

