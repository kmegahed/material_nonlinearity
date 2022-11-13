%clc; clear all;
lamda_target=360;dlamda_1=50;tol=1e-4;tol2=1e-4;ndofn=3;
for i=1:1
node=[0,0;
    0 288;
    240 288;
    720 288;
    720 0;
    240 0];
% (N1 N2) E,I, A, Fy, Zx
ele=[1 2 29000 248 13.3 36 54.9;
    6 3 29000 248 13.3 36 54.9;
    5 4 29000 248 13.3 36 54.9;
    2 3 29000 2850 24.8 36 244;
    3 4 29000 2850 24.8 36 244];
%node dof
bc=[1 1;1 2;1 3;
    5 1;5 2;5 3;
    6 1;6 2;6 3];
%node dof load
load=[2 1 0.1; 3 2 -1; 4 2 -0.5];
end
for i=1:1
    nel=size(ele,1);nnode=size(node,1);ndof=nnode*ndofn;node_updated=node;
    fi=zeros(ndofn*2,nel);Fi=fi;Fi_history=fi;up=fi;uT=fi;uT_history=fi;%(uT Total)
    du=zeros(ndof,1);u=du;u_history=du;pref=du;
    n=0;zz=1;lamda=0;hinge=[];KK=zeros(ndof);KE=zeros(6,nel);
    for j=1:size(load,1);pref(load(j,1)*ndofn-ndofn+load(j,2))=load(j,3);end
    const=[];for j=1:size(bc,1);const=[const;bc(j,1)*ndofn-ndofn+bc(j,2)];end; free=setdiff(1:ndof,const);
end
while lamda<=lamda_target
    if lamda==lamda_target;disp('max applied load');break;end
    if lamda+dlamda_1>lamda_target;dlamda=lamda_target-lamda;else; dlamda=dlamda_1;end
    Kt=KK;%zeros
    for iel=1:nel
        ele_nodes=ele(iel,1:2);nody=node_updated(ele_nodes,:);
        eldofs=         ndofn*(ele_nodes(1)-1)+1:ndofn*ele_nodes(1);
        eldofs=[eldofs  ndofn*(ele_nodes(2)-1)+1:ndofn*ele_nodes(2)];
        E=ele(iel,3); I=ele(iel,4); A=ele(iel,5); Py=A*ele(iel,6); Mp=ele(iel,6)*ele(iel,7);
        PM=[Fi(1,iel); Fi(3,iel); Fi(4,iel); Fi(6,iel)];%P1 M1 P2 M2
        [L,T]=trans(nody);ke_L=kee(E,L,A,I);kg_L=kgg(PM(3),L);
        [G,km_L]=kmm(ke_L,kg_L,PM,PM,Py,Mp,tol,0);%mid like previous PM
        kt_L=ke_L+kg_L+km_L;kt_G= T'*kt_L*T;Kt(eldofs,eldofs)=Kt(eldofs,eldofs)+kt_G;
    end
    keff=Kt(free,free);
    [V,D]=eig(keff); 
    if find(diag(D<=0));disp('illyy condition');break;
    elseif rcond(decomposition(keff))<eps;disp('illyy condition');break;
    end
    du_free=keff\(pref(free));
    tau_repeat_corrector=true;
    while tau_repeat_corrector
        tau_repeat=true;
        while tau_repeat
            du(free)=dlamda*du_free;
            node_temp=add_disp(node_updated,du,nnode,ndofn);
            dfint=KE;tau_array=[];
            for iel=1:nel
                ele_nodes=ele(iel,1:2);nody=node_updated(ele_nodes,:);
                eldofs=         ndofn*(ele_nodes(1)-1)+1:ndofn*ele_nodes(1);
                eldofs=[eldofs  ndofn*(ele_nodes(2)-1)+1:ndofn*ele_nodes(2)];
                E=ele(iel,3); I=ele(iel,4); A=ele(iel,5); Py=A*ele(iel,6); Mp=ele(iel,6)*ele(iel,7);
                PM=[Fi(1,iel); Fi(3,iel); Fi(4,iel); Fi(6,iel)];%P1 M1 P2 M2
                [L,T]=trans(nody);ke_L=kee(E,L,A,I);kg_L=kgg(PM(3),L);
                [G,km_L]=kmm(ke_L,kg_L,PM,PM,Py,Mp,tol,0);%mid like previous PM
                kt_L=ke_L+kg_L+km_L;
                du_e=T*du(eldofs);
                dfint(:,iel)=kt_L*du_e;fint(:,iel)=Fi(:,iel)+dfint(:,iel);
                [L,T_temp]=trans(node_temp(ele_nodes,:));fint(:,iel)=T_temp*T'*fint(:,iel);
                [tau_array,ph1_prev,ph2_prev]=tauy(tau_array,fint(:,iel),Fi(:,iel),Mp,Py,tol);
                if abs(ph1_prev-1)<=tol || abs(ph2_prev-1)<=tol
                    plastic_sign=(((G'*(ke_L+kg_L)*G)^-1)*G'*(ke_L+kg_L))*du_e;
                    if find(plastic_sign<0);disp('attemped unloading');return;end
                end
            end%element
            if find(tau_array)
                tt=min(tau_array);dlamda=tt*dlamda;
            else
                dfint_predictor=0.5*dfint;du_predictor=0.5*du;tau_repeat=false;
            end
        end %while tau_repeat
        node_predictor=add_disp(node_updated,du_predictor,nnode,ndofn);
        fint_predictor=Fi+dfint_predictor;
        for iel=1:nel
            ele_nodes=ele(iel,1:2);
            eldofs=         ndofn*(ele_nodes(1)-1)+1:ndofn*ele_nodes(1);
            eldofs=[eldofs  ndofn*(ele_nodes(2)-1)+1:ndofn*ele_nodes(2)];
            E=ele(iel,3); I=ele(iel,4); A=ele(iel,5); Py=A*ele(iel,6); Mp=ele(iel,6)*ele(iel,7);
            [L,T]=trans(node_updated(ele_nodes,:));
            [L,T2]=trans(node_predictor(ele_nodes,:));
            fint_predictor(:,iel)=T2*T'*fint_predictor(:,iel);
            p1=Fi(1,iel)/Py;p2=Fi(4,iel)/Py;m1=Fi(3,iel)/Mp;m2=Fi(6,iel)/Mp;
            ph1_prev=phii(p1,m1,0);ph2_prev=phii(p2,m2,0);
            p1_T=fint_predictor(1,iel)/Py;p2_T=fint_predictor(4,iel)/Py;m1_T=fint_predictor(3,iel)/Mp;m2_T=fint_predictor(6,iel)/Mp;
            ph1_T=phii(p1_T,m1_T,0);ph2_T=phii(p2_T,m2_T,0);ph1_tol=phii(p1_T,m1_T,0.01);ph2_tol=phii(p2_T,m2_T,0.01);
            if abs(ph1_prev-1)<=tol&&ph1_tol<=(1+tol)&&   (ph1_T>(1+tol) || ph1_T<(1-tol) )
                if ph1_T>(1+tol); flag=1;%disp('drifted point above yield surface.')
                elseif ph1_T<(1-tol); flag=0;%disp('drifted point below yield surface.')
                end
                [correc]=radial_projection(p1_T,m1_T,tol2,flag);
                if flag==1;fint_predictor([1 3],iel)=correc*fint_predictor([1 3],iel);
                else fint_predictor([1 3],iel)=fint_predictor([1 3],iel)+[correc(1)*Py;correc(2)*Mp];
                end
            end
            if abs(ph2_prev-1)<=tol&&ph2_tol<=(1+tol)&&   (ph2_T>(1+tol) || ph2_T<(1-tol) )
                if ph2_T>(1+tol); flag=1;%disp('drifted point above yield surface.')
                elseif ph2_T<(1-tol); flag=0;%disp('drifted point below yield surface.')
                end
                [correc]=radial_projection(p2_T,m2_T,tol2,flag);
                if flag==1;fint_predictor([4 6],iel)=correc*fint_predictor([4 6],iel);
                else fint_predictor([4 6],iel)=fint_predictor([4 6],iel)+[correc(1)*Py;correc(2)*Mp];
                end
            end
%             if abs(ph1_prev-1)<=tol&&ph1_tol<=(1+tol)&&   (ph1_T>(1+tol) || ph1_T<(1-tol) )
%                 [msol,psol,alphasol]=normal_projection(m1_T,p1_T,Py,Mp);
%                 fint_predictor([1 3],iel)=[psol*Py;msol*Mp];end
%             if abs(ph2_prev-1)<=tol&&ph2_tol<=(1+tol)&&   (ph2_T>(1+tol) || ph2_T<(1-tol) )
%                 [msol,psol,alphasol]=normal_projection(m2_T,p2_T,Py,Mp);
%                 fint_predictor([4 6],iel)=[psol*Py;msol*Mp];end
        end%element
        kt_corrector=KK;
        for iel=1:nel
            ele_nodes=ele(iel,1:2);nody=node_predictor(ele_nodes,:);
            eldofs=         ndofn*(ele_nodes(1)-1)+1:ndofn*ele_nodes(1);
            eldofs=[eldofs  ndofn*(ele_nodes(2)-1)+1:ndofn*ele_nodes(2)];
            E=ele(iel,3); I=ele(iel,4); A=ele(iel,5); Py=A*ele(iel,6); Mp=ele(iel,6)*ele(iel,7);
            PM=[Fi(1,iel); Fi(3,iel); Fi(4,iel); Fi(6,iel)];%P1 M1 P2 M2
            PM_P=[fint_predictor(1,iel); fint_predictor(3,iel); fint_predictor(4,iel); fint_predictor(6,iel)];
            [L,T]=trans(nody);ke_L=kee(E,L,A,I);kg_L=kgg(PM_P(3),L);
            [G,km_L]=kmm(ke_L,kg_L,PM,PM_P,Py,Mp,tol,0);%previous - mid (predictor) && control=0
            kt_L=ke_L+kg_L+km_L;kt_G= T'*kt_L*T;kt_corrector(eldofs,eldofs)=kt_corrector(eldofs,eldofs)+kt_G;
        end
        keff=kt_corrector(free,free);
        [V,D]=eig(keff);
        if find(diag(D<=0));disp('illyy condition');break;
        elseif rcond(decomposition(keff))<eps;disp('illyy condition');break;
        end
        du(free)=keff\(dlamda*pref(free));
        node_temp=add_disp(node_updated,du,nnode,ndofn);dfint=KE;tau_array=[];dup=KE;du_total=KE;
        for iel=1:nel
            ele_nodes=ele(iel,1:2);nody=node_predictor(ele_nodes,:);%predictor
            eldofs=         ndofn*(ele_nodes(1)-1)+1:ndofn*ele_nodes(1);
            eldofs=[eldofs  ndofn*(ele_nodes(2)-1)+1:ndofn*ele_nodes(2)];
            E=ele(iel,3); I=ele(iel,4); A=ele(iel,5); Py=A*ele(iel,6); Mp=ele(iel,6)*ele(iel,7);
            PM=[Fi(1,iel); Fi(3,iel); Fi(4,iel); Fi(6,iel)];%P1 M1 P2 M2
            PM_P=[fint_predictor(1,iel); fint_predictor(3,iel); fint_predictor(4,iel); fint_predictor(6,iel)];%new
            [L,T]=trans(nody);ke_L=kee(E,L,A,I);kg_L=kgg(PM_P(3),L);
            [G,km_L]=kmm(ke_L,kg_L,PM,PM_P,Py,Mp,tol,0);%same as last one
            kt_L=ke_L+kg_L+km_L;
            du_e=T*du(eldofs);
            dfint(:,iel)=kt_L*du_e;
            [L,T3]=trans(node_updated(ele_nodes,:));
            dfint(:,iel)=T3*T'*dfint(:,iel);
            fint(:,iel)=Fi(:,iel)+dfint(:,iel);
            [L,T2]=trans(node_temp(ele_nodes,:));
            fint(:,iel)=T2*T3'*fint(:,iel);
            [tau_array,ph1_prev,ph2_prev]=tauy(tau_array,fint(:,iel),Fi(:,iel),Mp,Py,tol);
            if abs(ph1_prev-1)<=tol || abs(ph2_prev-1)<=tol
                plastic_sign=(((G'*(ke_L+kg_L)*G)^-1)*G'*(ke_L+kg_L))*du_e;
                if find(plastic_sign<0);disp('attemped unloading');return;end
            end
        end%element
        if find(tau_array)
            tt=min(tau_array);dlamda=tt*dlamda;
        else
            for iel=1:nel
                ele_nodes=ele(iel,1:2);
                eldofs=         ndofn*(ele_nodes(1)-1)+1:ndofn*ele_nodes(1);
                eldofs=[eldofs  ndofn*(ele_nodes(2)-1)+1:ndofn*ele_nodes(2)];
                E=ele(iel,3); I=ele(iel,4); A=ele(iel,5); Py=A*ele(iel,6); Mp=ele(iel,6)*ele(iel,7);
%%%%%%%%%%%%%%%%%%%%%%
                [L,T3]=trans(node_updated(ele_nodes,:));
                [L,T]=trans(node_predictor(ele_nodes,:));
                [L,T2]=trans(node_temp(ele_nodes,:));
                dfint(:,iel)=T3*T'*dfint(:,iel);
                fint(:,iel)=Fi(:,iel)+dfint(:,iel);
                fint(:,iel)=T2*T3'*fint(:,iel);
%%%%%%%%%%%%%%%%%%%%%%%
                p1=Fi(1,iel)/Py;p2=Fi(4,iel)/Py;m1=Fi(3,iel)/Mp;m2=Fi(6,iel)/Mp;
                ph1_prev=phii(p1,m1,0);ph2_prev=phii(p2,m2,0);
                p1_T=fint(1,iel)/Py;p2_T=fint(4,iel)/Py;m1_T=fint(3,iel)/Mp;m2_T=fint(6,iel)/Mp;
                ph1_T=phii(p1_T,m1_T,0);ph2_T=phii(p2_T,m2_T,0);ph1_tol=phii(p1_T,m1_T,0.01);ph2_tol=phii(p2_T,m2_T,0.01);
                if abs(ph1_prev-1)<=tol&&ph1_tol<=(1+tol)&&   (ph1_T>(1+tol) || ph1_T<(1-tol) )
                    if ph1_T>(1+tol); flag=1;%disp('drifted point above yield surface.')
                    elseif ph1_T<(1-tol); flag=0;%disp('drifted point below yield surface.')
                    end
                    [correc]=radial_projection(p1_T,m1_T,tol2,flag);
                    if flag==1;fint([1 3],iel)=correc*fint([1 3],iel);
                    else fint([1 3],iel)=fint([1 3],iel)+[correc(1)*Py;correc(2)*Mp];
                    end
                end
                if abs(ph2_prev-1)<=tol&&ph2_tol<=(1+tol)&&   (ph2_T>(1+tol) || ph2_T<(1-tol) )
                    if ph2_T>(1+tol); flag=1;%disp('drifted point above yield surface.')
                    elseif ph2_T<(1-tol); flag=0;%disp('drifted point below yield surface.')
                    end
                    [correc]=radial_projection(p2_T,m2_T,tol2,flag);
                    if flag==1;fint([4 6],iel)=correc*fint([4 6],iel);
                    else fint([4 6],iel)=fint([4 6],iel)+[correc(1)*Py;correc(2)*Mp];
                    end
                end
%                 if abs(ph1_prev-1)<=tol&&ph1_tol<=(1+tol)&&   (ph1_T>(1+tol) || ph1_T<(1-tol) )
%                     [msol,psol,alphasol]=normal_projection(m1_T,p1_T,Py,Mp);
%                     fint([1 3],iel)=[psol*Py;msol*Mp];end
%                 if abs(ph2_prev-1)<=tol&&ph2_tol<=(1+tol)&&   (ph2_T>(1+tol) || ph2_T<(1-tol) )
%                     [msol,psol,alphasol]=normal_projection(m2_T,p2_T,Py,Mp);
%                     fint([4 6],iel)=[psol*Py;msol*Mp];end
            end%element
            tau_repeat_corrector=false;
        end%find(tau)
    end %while tau_repeat_corrector
    node_updated=node_temp;
    Fi=fint;Fi_history=[Fi_history;Fi];u=u+du;up=up+dup;
    %up_history=[up_histroy;up];
    lamda=lamda+dlamda;zz=zz+1;
end



function [L,T]=trans(nody)
E1 =nody(2,:) - nody(1,:);L = norm(E1);E1 = E1/L;
E2 = [-E1(2) E1(1)];Qrot = [E1 ; E2];Qrot(3,3) = 1;
T = [Qrot zeros(3); zeros(3) Qrot];
end
function ke=kee(E,L,A,I)
ke =    [ A*E/L            0               0       -A*E/L       0             0;
    0          12*E*I/(L^3)    6*E*I/(L^2)   0    -12*E*I/(L^3)   6*E*I/(L^2);
    0           6*E*I/(L^2)     4*E*I/L      0     -6*E*I/(L^2)    2*E*I/L;
    -A*E/L            0             0         A*E/L       0              0;
    0           -12*E*I/(L^3)   -6*E*I/(L^2)  0    12*E*I/(L^3)    -6*E*I/(L^2);
    0            6*E*I/(L^2)     2*E*I/L    0     -6*E*I/(L^2)      4*E*I/L ];
end
function kg=kgg(p,L)
kg=p/L*[1  0        0       -1   0      0;
    0  6/5     L/10      0  -6/5    L/10;
    0  L/10 (2*L^2)/15   0  -L/10 -L^2/30;
    -1   0       0        1    0      0;
    0  -6/5   -L/10      0  6/5   -L/10;
    0  L/10  -L^2/30     0 -L/10  (2*L^2)/15];
end
function ph=phii(p,m,t)
ph=(p/(1+t))^2+(m/(1+t))^2+3.5*(p/(1+t))^2*(m/(1+t))^2;
end
function [dph_p,dph_m]=dphii(p,m,t,Py,Mp)
    dph_p=(1/(1+t)^2)*((2*p/Py)   +   (1/(1+t)^2)*(7*p  *m^2/Py));
    dph_m=(1/(1+t)^2)*((2*m/Mp)   +   (1/(1+t)^2)*(7*p^2*m  /Mp));
end
function tau_r=falsi_midstep_surface(p,m,tol)
converge=false;phi_trial_mid=phii(p,m,0);
if phi_trial_mid>(1+tol);    tau_a=0;tau_b= 0.01;
elseif phi_trial_mid<(1-tol);tau_a=0;tau_b=-0.01;
else                         tau_r=0;converge=true;
end

while ~converge
    phi_a=phii(p,m,tau_a);phi_b=phii(p,m,tau_b);    
    tau_r= tau_b - ((phi_b -1)*((tau_a - tau_b)/(phi_a - phi_b)));
    phi_r=phii(p,m,tau_r);    
    if abs(phi_r -1) <= tol;    converge=true;break;end
    if (phi_r -1)*(phi_b -1)>0;     tau_b=tau_r;
    elseif (phi_r-1)*(phi_a -1)>0;  tau_a=tau_r;
    end
end
end
function [G,km]=kmm(ke,kg,fp,fm,Py,Mp,tol,control)%fp previous fm mean mid
p1_prev=fp(1)/Py;m1_prev=fp(2)/Mp;  p2_prev=fp(3)/Py;m2_prev=fp(4)/Mp;
p1_mid=fm(1)/Py;m1_mid=fm(2)/Mp;    p2_mid=fm(3)/Py;m2_mid=fm(4)/Mp;
phi_1_prev=phii(p1_prev,m1_prev,0);phi_2_prev=phii(p2_prev,m2_prev,0);
if control==1
    if abs(phi_1_prev -1)<=tol
        [tau_r]=falsi_midstep_surface(p1_mid,m1_mid,tol);
        [dphip1,dphim1]=dphii(p1_mid,m1_mid,tau_r,Py,Mp);
    end
    if abs(phi_2_prev -1)<=tol
        [tau_r]=falsi_midstep_surface(p2_mid,m2_mid,tol);
        [dphip2,dphim2]=dphii(p2_mid,m2_mid,tau_r,Py,Mp);
    end
else
    if abs(phi_1_prev -1)<=tol
        [dphip1,dphim1]=dphii(p1_mid,m1_mid,0,Py,Mp);
    end
    if abs(phi_2_prev -1)<=tol
        [dphip2,dphim2]=dphii(p2_mid,m2_mid,0,Py,Mp);
    end
end
if phi_1_prev<(1-tol) && phi_2_prev<(1-tol)             G=0;km=zeros(6);
else
    if abs(phi_1_prev-1)<=tol && phi_2_prev<(1-tol)
        G=[dphip1;0;dphim1;0;0;0];
    elseif phi_1_prev<(1-tol) && abs(phi_2_prev-1)<=tol
        G=[0;0;0;dphip2;0;dphim2;];
    elseif abs(phi_1_prev-1)<=tol && abs(phi_2_prev-1)<=tol
        G=[dphip1 0;
            0     0;
            dphim1 0;
            0 dphip2;
            0     0;
            0 dphim2];
    end
    km=-(ke+kg)*G*((G'*(ke+kg)*G)^-1)*G'*(ke+kg);
end
end
function [msol,psol,alphasol]=normal_projection(m1,p1,Py,Mp)

%Welcome to drift corrector routine. This script is written to return trial
%fint values (that have been scaled to remain within expanded surface) back
%to original yield surface. Corrected fint trial values are returned via
%gradient to the original surface as opposed to radial return used
%previously. Matlab's vpasolve is invoked to solve three nonlinear
%equations for m, p, alpha parameters. 

syms m p alp real
eq1=m*(2*alp + 7*alp*p^2 + Mp)== m1*Mp;
eq2=p*(2*alp + 7*alp*m^2 + Py) == p1*Py;
eq3=m^2 + p^2 + 3.5*p^2*m^2==1;
%[sol1,sol2,sol3]=vpasolve(eq1,eq2,eq3,[m,p,alp]);
%S=vpasolve(eq1,eq2,eq3,m,[0 1],p,[0 1]);
vars=[m p alp];range=[0 sign(m1)*1;0 sign(p1)*1;NaN NaN];
S=vpasolve(eq1,eq2,eq3,vars,range);
msol=S.m;psol=S.p;alphasol=S.alp;
end
function node_temp=add_disp(node_updated,du,nnode,ndofn)
       for i=1:nnode
                node_temp(i,1)=node_updated(i,1)+du(i*ndofn-2);node_temp(i,2)=node_updated(i,2)+du(i*ndofn-1);
       end
end
function tau=elastic_return(PM_prev,PM_trial,Py,Mp,tol,casy)
p=PM_prev(1)/Py;m=PM_prev(2)/Mp;dp=(PM_trial(1)/Py)-p;dm=(PM_trial(2)/Mp)-m;
tau_a=0;tau_b=1;        converge=false;j=1;
while ~converge %|| j==50
    phi_a=phii(p+tau_a*dp,m+tau_a*dm,0.01*(casy-1));phi_b=phii(p+tau_b*dp,m+tau_b*dm,0.01*(casy-1));
    tau_r=tau_b - ((phi_b -1)*(tau_a-tau_b))/(phi_a - phi_b);
    phi_r=phii(p+tau_r*dp,m+tau_r*dm,0.01*(casy-1));
    if abs(phi_r-1)<=tol;        converge=true;        break;    end
    if (phi_r -1)*(phi_b-1)>0;tau_b=tau_r;elseif (phi_r -1)*(phi_a -1)>0; tau_a=tau_r; end
    j=j+1;
end
tau=tau_r;
end
function [tau_array,ph1_prev,ph2_prev]=tauy(tau_array,fint,Fi,Mp,Py,tol)
p1=Fi(1)/Py;p2=Fi(4)/Py;m1=Fi(3)/Mp;m2=Fi(6)/Mp;
ph1_prev=phii(p1,m1,0);ph2_prev=phii(p2,m2,0);
p1T=fint(1)/Py;p2T=fint(4)/Py;m1T=fint(3)/Mp;m2T=fint(6)/Mp;
ph1_T=phii(p1T,m1T,0);ph2_T=phii(p2T,m2T,0);
ph1_tol=phii(p1T,m1T,0.01);ph2_tol=phii(p2T,m2T,0.01);
ph1_alarm=phii(p1T,m1T,-0.01);ph2_alarm=phii(p2T,m2T,-0.01);
if (ph1_alarm<=(1+tol)&&abs(ph1_prev-1)<=tol)||((ph2_alarm<=(1+tol))&&abs(ph2_prev-1)<=tol)
    %this means that trial pm point has gone under severe
    %degeneracy i.e. tangent vector to the original surface is
    %degenerate! analysis must terminate
    disp('Plastic flow vectors for some member ends might not be positioned along the surface properly. Internal force recovery might not be accurate.');
end
if ph1_T>1+tol&& ph1_prev<1-tol;casy=1;%elastic return
elseif ph1_tol>1+tol&&abs(ph1_prev-1)<=tol;casy=2; %drift control
else
    casy=0;
end
if casy~=0
    tau=elastic_return(Fi([1 3]),fint([1 3]),Py,Mp,tol,casy);
    tau_array=[tau_array tau];
end
if ph2_T>1+tol&& ph2_prev<1-tol;casy=1;
elseif ph2_tol>1+tol&&abs(ph2_prev-1)<=tol;casy=2; 
else
    casy=0;
end
if casy~=0
    tau=elastic_return(Fi([4 6]),fint([4 6]),Py,Mp,tol,casy);
    tau_array=[tau_array tau];
end
end
function [correc]=radial_projection(pvalue,mvalue,tol2,flag)
if flag==1
    p=0;m=0;dp=pvalue;dm=mvalue;tau_a=0;tau_b=1;
else
    p=pvalue;m=mvalue;dp=sqrt(pvalue^2 + mvalue^2)*pvalue;
    dm=sqrt(pvalue^2 + mvalue^2)*mvalue;
    tau_a=0;tau_b=1;
end
%tol=0.01; %tolerance parameter for convergence to phi==1
converge=false;j=1;
while ~converge %|| j==50
    phi_a=phii(p+tau_a*dp,m+tau_a*dm,0);phi_b=phii(p+tau_b*dp,m+tau_b*dm,0);
    tau_r=tau_b - ((phi_b -1)*(tau_a-tau_b))/(phi_a - phi_b);
    phi_r=phii(p+tau_r*dp,m+tau_r*dm,0);
    if (abs(phi_r-1))<=tol2; converge=true;break;end
    if (phi_r -1)*(phi_b-1)>0;          tau_b=tau_r;
    elseif (phi_r -1)*(phi_a -1)>0;     tau_a=tau_r;
    end
    j=j+1;
end
if flag==1
    correc=tau_r;
else
    correc=tau_r*[dp;dm];
end
end
