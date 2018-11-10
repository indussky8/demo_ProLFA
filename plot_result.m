
% plot(X(:,1),X(:,2),'p','LineWidth',1,'color','b','markersize',3)
% hold on
% Pro = X'*Z;  %d*c cprototypes
% plot(Pro(1,:),Pro(2,:),'o','LineWidth',1,'color','r','markersize',7)
% legend('descriptors','prototypes','FontSize',24,'FontName','Times New Roman')
% %title('Source set','FontSize',24,'FontName','Times New Roman')
% set(gca,'FontSize',10)

% %  
% X = X';
% plot(X(1,:),X(2,:),'p','LineWidth',1,'color','b','markersize',3)
% hold on
% plot(b(:,1),b(:,2),'o','LineWidth',1,'color','r','markersize',7)
% legend('descriptors','prototypes','FontSize',24,'FontName','Times New Roman')
% %title('Source set','FontSize',24,'FontName','Times New Roman')
% set(gca,'FontSize',10)

% 

plot(X(1:20,1),X(1:20,2),'p','LineWidth',1,'color','b','markersize',3)
% hold on
% % plot(X(21:40,1),X(21:40,2),'p','LineWidth',1,'color','b','markersize',3)
% % hold on
% % plot(X(41:60,1),X(41:60,2),'p','LineWidth',1,'color','b','markersize',3)
% % hold on
% % plot(X(61:80,1),X(61:80,2),'p','LineWidth',1,'color','b','markersize',3)
% % hold on
% % plot(X(81:104,1),X(81:104,2),'p','LineWidth',1,'color','b','markersize',3)
% % hold on
% % plot(X(105:120,1),X(105:120,2),'^','LineWidth',1,'color','c','markersize',3)
% % hold on
% % plot(X(121:140,1),X(121:140,2),'^','LineWidth',1,'color','c','markersize',3)
% % hold on
% % plot(X(141:160,1),X(141:160,2),'^','LineWidth',1,'color','c','markersize',3)
% % hold on
% % plot(X(161:180,1),X(161:180,2),'^','LineWidth',1,'color','c','markersize',3)
% % hold on
% % plot(X(181:200,1),X(181:200,2),'^','LineWidth',1,'color','c','markersize',3)
% % hold on
% 
plot(X(1:104,1),X(1:104,2),'p','LineWidth',1,'color','g','markersize',3)
hold on
plot(X(105:200,1),X(105:200,2),'^','LineWidth',1,'color','c','markersize',3)
hold on


Prototype = X'*Z_200;
plot(Prototype(1,:),Prototype(2,:),'+','LineWidth',1,'color','k','markersize',7)

hold on
Pro = G*X*X'*Z;  %d*c cprototypes
% Pro = G*X*b';  %d*c cprototypes
plot(Pro(1:5,1),Pro(1:5,2),'s','LineWidth',1,'color','k','markersize',7)
hold on
plot(Pro(6:10,1),Pro(6:10,2),'o','LineWidth',1,'color','k','markersize',7)
legend('local descriptor in the 1st class','local descriptor in the 2nd class','the selected 2 prototypes','global representation of each sample in the 1st class','global representation of each sample in the 2nd class','FontSize',24,'FontName','Times New Roman')
%title('Source set','FontSize',24,'FontName','Times New Roman')


%   plot(X(45,:),Prototype(2,:),'+','LineWidth',1,'color','k','markersize',7)
















set(gca,'FontSize',10)
