function Processing

global timeB B



min(diff(timeB))
max(diff(timeB))

%half a second PVI
dt=min(diff(timeB));    %minimum time step
Ntau=floor(0.2./dt);    %# of point to compute PVI

PVIx=B{1}((Ntau+1):end)-B{1}(1:end-Ntau);
PVIy=B{2}((Ntau+1):end)-B{2}(1:end-Ntau);
PVIz=B{3}((Ntau+1):end)-B{3}(1:end-Ntau);
PVI=sqrt(PVIx.^2+PVIy.^2+PVIz.^2);
time_PVI=timeB(1:end-Ntau);


%min and max To values
tmin=min(timeB); tmax=max(timeB);
%time step
Jo={}; dT=20*60;
[I,Ja]=find(time_PVI'>=tmin+0*dT); [I,Jb]=find(time_PVI'<=tmin+1*dT); Jo{1}=intersect(Ja,Jb);
[I,Ja]=find(time_PVI'>=tmin+1*dT); [I,Jb]=find(time_PVI'<=tmin+2*dT); Jo{2}=intersect(Ja,Jb);
[I,Ja]=find(time_PVI'>=tmin+2*dT); [I,Jb]=find(time_PVI'<=tmin+3*dT); Jo{3}=intersect(Ja,Jb);
[I,Ja]=find(time_PVI'>=tmin+3*dT); [I,Jb]=find(time_PVI'<=tmin+4*dT); Jo{4}=intersect(Ja,Jb);
[I,Ja]=find(time_PVI'>=tmin+4*dT); [I,Jb]=find(time_PVI'<=tmin+5*dT); Jo{5}=intersect(Ja,Jb);
[I,Ja]=find(time_PVI'>=tmin+5*dT); [I,Jb]=find(time_PVI'<=tmax);     Jo{6}=intersect(Ja,Jb);



norm_PVIx=zeros(size(PVIx));
norm_PVIy=zeros(size(PVIy));
norm_PVIz=zeros(size(PVIz));

sigma_x=[]; sigma_y=[]; sigma_z=[];
for i=1:6
    range=Jo{i};
    time_sigma=mean(time_PVI(range));
    sigma_x=[sigma_x sqrt(mean(PVIx(range).^2))];
    sigma_y=[sigma_y sqrt(mean(PVIy(range).^2))];
    sigma_z=[sigma_z sqrt(mean(PVIz(range).^2))];
    
    norm_PVIx(range)=PVIx(range)./sigma_x(end);
    norm_PVIy(range)=PVIy(range)./sigma_y(end);
    norm_PVIz(range)=PVIz(range)./sigma_z(end);
end


norm_PVI=sqrt(norm_PVIx.^2+norm_PVIy.^2+norm_PVIz.^2);



figure('position',[100 100 1100 620]);
subplot(4,1,1);
hold on;
h=plot(time_PVI,norm_PVIx,'k-');
set(h,'linewidth',1);
set(gca,'fontsize',14);
set(gca,'xticklabel',[]);
xlim([tmin tmax]);
ylim([-10 10]);
ylabel('PVIx');
grid on;
title('PVI [0.2s]');

subplot(4,1,2);
hold on;
h=plot(time_PVI,norm_PVIy,'r-');
set(h,'linewidth',1);
set(gca,'fontsize',14);
set(gca,'xticklabel',[]);
xlim([tmin tmax]);
ylim([-10 10]);
ylabel('PVIy');
grid on;

subplot(4,1,3);
hold on;
h=plot(time_PVI,norm_PVIz,'b-');
set(h,'linewidth',1);
set(gca,'fontsize',14);
set(gca,'xticklabel',[]);
xlim([tmin tmax]);
ylim([-10 10]);
ylabel('PVIz');
grid on;

subplot(4,1,4);
hold on;
h=plot(time_PVI,norm_PVI,'k-');
set(h,'linewidth',0.5);
h=plot(time_PVI,5+zeros(size(time_PVI)),'r--');
set(h,'linewidth',2);
set(gca,'fontsize',14);
%set(gca,'yscale','log');
legend('PVI','PVI=5');
xlim([tmin tmax]);
ylim([0 10]);
xlabel('time');
ylabel('PVI');
grid on;


[I,num5]=find(norm_PVI'>=5);
[I,j5]=find(diff(num5)>1);
N5=length(j5);
interv={};
interv{1}=num5(1):num5(j5(1));
for i=1:(N5-1)
    interv{i+1}=num5(j5(i)+1):num5(j5(i+1));
    %a=interv{i+1}
end
interv{N5+1}=num5(j5(N5)+1):num5(end);
%a=interv{N5+1}


dir='plots 06-08\';
for i=1:(N5+1)
    proc=i/(N5+1)*100

    t_left=timeB(min(interv{i}));
    t_righ=timeB(max(interv{i}));
    to=(t_righ+t_left)/2;
    
    dt=max((t_righ-t_left)/2,2);
    [I,J1]=find(timeB'>=to-dt); [I,J2]=find(timeB'<=to+dt); Jo=intersect(J1,J2);
    
    %subselect magnetic fields
    tlmn=timeB(Jo);
    Bx=B{1}(Jo); By=B{2}(Jo); Bz=B{3}(Jo);
    %MVA coordinates
    [dirL,dirM,dirN,EigValue]=MVA(Bx,By,Bz);
    
    Bl=Bx*dirL(1)+By*dirL(2)+Bz*dirL(3);
    Bm=Bx*dirM(1)+By*dirM(2)+Bz*dirM(3);
    Bn=Bx*dirN(1)+By*dirN(2)+Bz*dirN(3);
    
    hf=figure('position',[50 50 800 600],'visible','off');
    subplot(2,1,1);
    hold on;
    h=plot(tlmn-to,Bl,'k-');
    set(h,'linewidth',2);
    h=plot(tlmn-to,Bm,'r-');
    set(h,'linewidth',2);
    h=plot(tlmn-to,Bn,'b-');
    set(h,'linewidth',2);
    set(gca,'fontsize',14);
    legend('B_{l}','B_{m}','B_{n}');
    set(gca,'xticklabel',[]);
    xlim([min(tlmn)-to max(tlmn)-to]);
    %xlabel('time [s]');
    ylabel('B [nT]');
    grid on;
    
    subplot(2,1,2);
    hold on;
    h=plot(tlmn-to,Bl-mean(Bl),'k-');
    set(h,'linewidth',2);
    h=plot(tlmn-to,Bm-mean(Bm),'r-');
    set(h,'linewidth',2);
    h=plot(tlmn-to,Bn-mean(Bn),'b-');
    set(h,'linewidth',2);
    set(gca,'fontsize',14);
    %set(gca,'xticklabel',[]);
    %legend('B_{l}-<B_{l}>','B_{m}-<B_{m}>','B_{n}-<B_{n}>');
    xlim([min(tlmn)-to max(tlmn)-to]);
    name=['t - ',num2str(to),' [s]'];
    xlabel(name);
    ylabel('B-<B_{\alpha}> [nT]');
    grid on;
    name=[dir,'i=',num2str(i),'.png'];
    saveas(hf,name);
    delete(hf);   
end



%{
figure;
hold on;
h=plot(timeB,B{1},'k-');
set(h,'linewidth',2);
h=plot(timeB,B{2},'r-');
set(h,'linewidth',2);
h=plot(timeB,B{3},'b-');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlabel('time');
ylabel('B [nT]');
%}


