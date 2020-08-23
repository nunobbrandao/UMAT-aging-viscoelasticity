clear all, close all,clc
set(0,'DefaultAxesFontname','tempos New Roman',...
'DefaultAxesFontsize',16,...
'DefaultTextFontSize', 18, ...
'DefaultTextFontWeight','bold',...
'DefaultLineLineWidth', 1.5, ...
'DefaultAxesFontname','Times New Roman',...
'DefaultTextFontname','Times New Roman')
#-----------------------------Import txt. files --------------------------------
%TIME(1); MISES(2); PRESSURE(3); COORDINATEXX_EXTERNAL_CASING(4); COORDINATEXX_INTERNAL_SALT(5)
filename = 'Node_data_P2001.dat';
delimiterIn = ' ';
headerlinesIn = 2;
Node = importdata(filename,delimiterIn,headerlinesIn);
Node_data_P2001=Node.data;
%NODE(1); COORDINATE(2); MISES(3+(i-1)*NUMBER_OF_FRAMES); PRESSURE(3+(i-1)*NUMBER_OF_FRAMES); (...)
filename = 'Path_data_P2001.dat';
delimiterIn = ' ';
headerlinesIn = 1;
Path = importdata(filename,delimiterIn,headerlinesIn);
Path_data_P2001=Path.data;
%TIME(1); MISES(2); PRESSURE(3); COORDINATEXX_EXTERNAL_CASING(4); COORDINATEXX_INTERNAL_SALT(5)
filename = 'Node_data_GL2007.dat';
delimiterIn = ' ';
headerlinesIn = 2;
Node = importdata(filename,delimiterIn,headerlinesIn);
Node_data_GL2007=Node.data;
%NODE(1); COORDINATE(2); MISES(3+(i-1)*NUMBER_OF_FRAMES); PRESSURE(3+(i-1)*NUMBER_OF_FRAMES); (...)
filename = 'Path_data_GL2007.dat';
delimiterIn = ' ';
headerlinesIn = 1;
Path = importdata(filename,delimiterIn,headerlinesIn);
Path_data_GL2007=Path.data;
%TIME(1); MISES(2); PRESSURE(3); COORDINATEXX_EXTERNAL_CASING(4); COORDINATEXX_INTERNAL_SALT(5)
filename = 'Node_data_R2015.dat';
delimiterIn = ' ';
headerlinesIn = 2;
Node = importdata(filename,delimiterIn,headerlinesIn);
Node_data_R2015=Node.data;
%NODE(1); COORDINATE(2); MISES(3+(i-1)*NUMBER_OF_FRAMES); PRESSURE(3+(i-1)*NUMBER_OF_FRAMES); (...)
filename = 'Path_data_R2015.dat';
delimiterIn = ' ';
headerlinesIn = 1;
Path = importdata(filename,delimiterIn,headerlinesIn);
Path_data_R2015=Path.data;
%----------------------------- INICIATE VARIAB: --------------------------------
pwc=8.9600E+04          %slurry hydrostatic pressure [kN/m2]
Path_data_FRAMES=4      %number of chosen frames in the Path_data
Path_data_QUANTITIES=6  %number of chosen quantities in the Path_data (mises, pressure, etc)
Path.colheaders         %show what and how many times I am plotting 
matrix={Node_data_P2001 Path_data_P2001 Node_data_GL2007 Path_data_GL2007 Node_data_R2015 Path_data_R2015};
leg={'P2001' 'GL2007' 'R2015'};
%
%---------------------PLOT: -------------------------
%
%PLOT1 - TIME VS ANNULUS FOR ONE NODE
%PLOT2 - TIME VS PRESSURE AND MISES FOR ONE NODE
%PLOT3 - COORDINATE VS PRESSURE AND MISES FOR # OF Path_data_frames 
%PLOT4 - COORDINATE VS STRESSES FOR # OF Path_data_frames 
%
n_analyses=size(matrix,2)/2
for j=1:n_analyses
  r=1;
  while (matrix{1+(j-1)*2}(r,3)~=pwc)
    r=r+1;%r represents the line when pwc is reached
  end
  last_line=size(matrix{1+(j-1)*2},1)
  figure(2); %plots Time vs Pressure and Mises 
    hold on                                                                            
    p(j)=plot(matrix{1+(j-1)*2}(r:last_line,1),matrix{1+(j-1)*2}(r:last_line,3),'k-','LineWidth',2,'Color', [j/(1+n_analyses) j/(1+n_analyses) j/(1+n_analyses)],'DisplayName',[leg{j} '-Press'] );
    m(j)=plot(matrix{1+(j-1)*2}(r:last_line,1),matrix{1+(j-1)*2}(r:last_line,2),'k-.','LineWidth',2,'Color',[j/(1+n_analyses) j/(1+n_analyses) j/(1+n_analyses)],'DisplayName',[leg{j} '-Mises']);
    for i=1:1:Path_data_FRAMES 
      figure(3);
        hold on
        cM(j,i)=plot(matrix{2+(j-1)*2}(2:100,2),matrix{2+(j-1)*2}(2:100,3+(i-1)*Path_data_QUANTITIES),'k-','LineWidth',2,'DisplayName', ['Mises' leg{j} '/' Path.colheaders{i+2}],'Color', [j/(1+n_analyses) j/(1+n_analyses) i/(1+n_analyses)]);
        cP(j,i)=plot(matrix{2+(j-1)*2}(2:100,2),matrix{2+(j-1)*2}(2:100,4+(i-1)*Path_data_QUANTITIES),'r-.','LineWidth',2,'DisplayName',['Press' leg{j} '/' Path.colheaders{i+2}],'Color', [j/(1+n_analyses) j/(1+n_analyses) i/(1+n_analyses)]);
      figure(4);   
        hold on
        s1(j,i)=plot(matrix{2+(j-1)*2}(2:100,2),matrix{2+(j-1)*2}(2:100,5+(i-1)*Path_data_QUANTITIES),'k-','LineWidth',2,'DisplayName',[leg{j} '/' Path.colheaders{i+2}],'Color', [j/(1+n_analyses) j/(1+n_analyses) i/(1+n_analyses)]);
      figure(5);   
        hold on      
        s2(j,i)=plot(matrix{2+(j-1)*2}(2:100,2),matrix{2+(j-1)*2}(2:100,6+(i-1)*Path_data_QUANTITIES),'r-','LineWidth',2,'DisplayName',[leg{j} '/' Path.colheaders{i+2}],'Color', [j/(1+n_analyses) j/(1+n_analyses) i/(1+n_analyses)]);
      figure(6);   
        hold on     
        s3(j,i)=plot(matrix{2+(j-1)*2}(2:100,2),matrix{2+(j-1)*2}(2:100,7+(i-1)*Path_data_QUANTITIES),'g-','LineWidth',2,'DisplayName',[leg{j} '/' Path.colheaders{i+2}],'Color', [j/(1+n_analyses) j/(1+n_analyses) i/(1+n_analyses)]);
      figure(7);   
        hold on      
        s4(j,i)=plot(matrix{2+(j-1)*2}(2:100,2),matrix{2+(j-1)*2}(2:100,8+(i-1)*Path_data_QUANTITIES),'b-','LineWidth',2,'DisplayName',[leg{j} '/' Path.colheaders{i+2}],'Color', [j/(1+n_analyses) j/(1+n_analyses) i/(1+n_analyses)]);
    end
end
%---------------------FORMAT PLOTS: -------------------------
figure(2)
    legend([p m],'Location','eastoutside','orientation','vertical');
    legend show 
    limy=[0:1e4:10e4]; 
    limx=[31.105:2:31.105+28];  
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Time [days]')
    ylabel('stress [MPa]')
    set(gca,'XTick',limx) 
    set(gca,'YTick',limy) 
    set(gca,'YTickLabel',sprintf('%0.0f|',limy/1000))
    set(gca,'XTickLabel',sprintf('%0.0f|',limx-1.105))
    box on 
    print -dpdf -color figure1.pdf
    limx=[31.105:.1:31.105+1];
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Time [minutes]')
    set(gca,'XTick',limx) 
    set(gca,'XTickLabel',sprintf('%0.0f|',(limx-31.105)*24*60)) 
    print -dpdf -color figure2.pdf 
figure(3)
    legend([cM cP],'Location','eastoutside','orientation','vertical');
    legend show   
    limy=[0:5e4:30e4]; 
    limx=[.15:0.03:.50];  
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Coordi. [cm]')
    ylabel('stress [MPa]')
    set(gca,'XTick',limx) 
    set(gca,'YTick',limy) 
    set(gca,'YTickLabel',sprintf('%0.0f|',limy/1000))
    set(gca,'XTickLabel',sprintf('%0.0f|',limx*100))
    box on
    print -dpdf -color figure3.pdf 
figure(4)
    legend(s1,'Location','eastoutside','orientation','vertical');
    legend show   
    limy=[-100E3:5e3:-5e4];  
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Coordi. [cm]')
    ylabel('stress [MPa]')
    set(gca,'XTick',limx) 
    set(gca,'YTick',limy) 
    set(gca,'YTickLabel',sprintf('%0.0f|',limy/1000))
    set(gca,'XTickLabel',sprintf('%0.0f|',limx*100))
    grid on
    box on
    title('S11')
    print -dpdf -color figure4.pdf 
figure(5)
    legend(s2,'Location','eastoutside','orientation','vertical');
    legend show   
    limy=[-400E3:5e4:-5e4];  
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Coordi. [cm]')
    ylabel('stress [MPa]')
    set(gca,'XTick',limx) 
    set(gca,'YTick',limy) 
    set(gca,'YTickLabel',sprintf('%0.0f|',limy/1000))
    set(gca,'XTickLabel',sprintf('%0.0f|',limx*100))
    grid on
    box on
    title('S22')
    print -dpdf -color figure5.pdf 
figure(6)
    legend(s3,'Location','eastoutside','orientation','vertical');
    legend show   
    limy=[-130E3:1e4:-70e3];  
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Coordi. [cm]')
    ylabel('stress [MPa]')
    set(gca,'XTick',limx) 
    set(gca,'YTick',limy) 
    set(gca,'YTickLabel',sprintf('%0.0f|',limy/1000))
    set(gca,'XTickLabel',sprintf('%0.0f|',limx*100))
    grid on
    box on
    title('S33')
    print -dpdf -color figure6.pdf 
figure(7)
    legend(s4,'Location','eastoutside','orientation','vertical');
    legend show   
    limy=[-0.4E3:0.4e3:2e3];  
    axis([min(limx) max(limx) min(limy) max(limy)])
    xlabel('Coordi. [cm]')
    ylabel('stress [MPa]')
    set(gca,'XTick',limx) 
    set(gca,'YTick',limy) 
    set(gca,'YTickLabel',sprintf('%0.2f|',limy/1000))
    set(gca,'XTickLabel',sprintf('%0.0f|',limx*100))
    grid on
    box on
    title('S12')
    print -dpdf -color figure7.pdf 
    