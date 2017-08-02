%% define global variable and load the stagnation point data

global vid h h1 h2 h9 S_N S_S S_E S_W options dim theta
options=optimset('TolX',1e-10,'TolFun',1e-10,'Algorithm','levenberg-marquardt','Display','off'); 

%% clear video and workspace data

clc;
close all;
clearvars -except center dim theta
imaqreset

try
closepreview();
stop(vid);
flushdata(vid);
delete(vid);
clear vid

imaqreset
catch ME
end
warning('off','all');

%% Define other operation variables

i=1;
thresh=0.2;
frames=zeros(1200,1600);
IFT_track = 0;

%% create a  a video out

vid=videoinput('qimaging',1);

set(vid,'FramesPerTrigger', inf);
triggerconfig(vid,'manual');
vid.Timeout = 30;

src = getselectedsource(vid);
set(src,'Exposure',0.100);

%% determine the device center

try
    if isempty(center)
        [center,dim,theta]= device_center(vid, thresh);
    end
catch
    [center,dim,theta]= device_center(vid, thresh);
end 

%% define field of view 

options=optimset('TolX',1e-10,'TolFun',1e-10,'Algorithm','levenberg-marquardt','Display','off');

solxy=[0.5,0.5]*dim;

global imref 
imref=(15*getsnapshot(vid));
low_in = 799/65535;
high_in = 6181/65535;
imref = imadjust(imref,[low_in high_in]);

%% Initial inlet and outlet pressures

P_N =1;
P_W =0.20;
P_S =1;
P_E =0.20;

V1=P_S;
V2=P_W;
V3=P_N;
V4=P_E;

%% create a 'NI' data session

s = daq.createSession('ni');
s.Rate=10000;

s.addAnalogOutputChannel('Dev8',0,'Voltage'); %%% Input 2 S %%%
s.addAnalogOutputChannel('Dev6',0,'Voltage'); %%% Output 1 E %%%
s.addAnalogOutputChannel('Dev8',1,'Voltage'); %%% Input 1 N %%%
s.addAnalogOutputChannel('Dev7',0,'Voltage'); %%% Output 2 W %%%

s.outputSingleScan([V1,V2,V3,V4]);

%% Define variables to save image, voltage and elapsed time history

im=uint16(zeros(1200,1600,1));
Voltage_data=zeros(1,13);
data=zeros(1,13);

%% Interactive consol
        
        h=figure(10); set(h,'units','normalized','outerposition',[0 0 1 1]);
        
        h1 = axes('Parent',h,'Units','normalized','Position',[0.5 0.22 0.45 .72]);axis off;

        h2= uicontrol('Parent',h,'units','normalized','Style', 'popup','String', 'Manual Control |IFT|Coalescence|Xs play|Xs Control|exit','Position', [0.55 0.2 0.1 0.02],'value',1);
        Mode_input=uicontrol('Parent',h,'units','normalized','Style', 'popup','Position', [0.7 0.1 0.1 0.02],'string','Proportional control|Stagnation point control');

        Exposure_text=uicontrol('Parent',h,'units','normalized','Style', 'text','Position', [0.85 0.2 0.04 0.02],'string','Exposure:','FontSize',10);
        Exposure=uicontrol('Parent',h,'units','normalized','Style', 'edit','Position', [0.9 0.2 0.025 0.02],'string','0.1');
        Toggle_dropdown_box=uicontrol('Parent',h,'units','normalized','Style', 'popup','Position', [0.55 0.1 0.1 0.02],'string','Click mode OFF|Click mode ON');
        Gain_text=uicontrol('Parent',h,'units','normalized','Style', 'text','Position', [0.85 0.15 0.04 0.02],'string','Gain:','FontSize',10);
        Gain=uicontrol('Parent',h,'units','normalized','Style', 'edit','Position', [0.9 0.15 0.025 0.02],'string','108');
        
        symbols_on=uicontrol('Parent',h,'units','normalized','Style', 'checkbox','Position', [0.7 0.15 0.1 0.02],'string','Symbols On','value',1);
      
        ift_image_saving_on=uicontrol('Parent',h,'units','normalized','Style', 'checkbox','Position', [0.7 0.20 0.1 0.02],'string','IFT Saving On','value',1);
      

        P_MIN_text=uicontrol('Parent',h,'units','normalized','Style', 'text','Position', [0.85 0.1 0.04 0.02],'string','P_MIN:','FontSize',10);
        P_MIN=uicontrol('Parent',h,'units','normalized','Style', 'edit','Position', [0.9 0.1 0.025 0.02],'string','0.5');
        
        GainIFT_text=uicontrol('Parent',h,'units','normalized','Style', 'text','Position', [0.85 0.05 0.04 0.02],'string','GainIFT:','FontSize',10);
        GainIFT=uicontrol('Parent',h,'units','normalized','Style', 'edit','Position', [0.9 0.05 0.025 0.02],'string','0.003');

        MC_panel=uipanel('Parent',h,'Title','Manual Control','background','white','FontSize',10,'Position',[0.03 .55 .4 .35]);

        S_N=uicontrol('Parent',MC_panel,'units','normalized','position',[0.05,0.75,0.9,0.05],'style','slider','Min',0.0,'Max',5.0,'value',P_N);
        Max_N=uicontrol('Parent',MC_panel,'units','normalized','position',[0.92,0.7,0.05,0.04],'style','text','string','5.0');
        Min_N=uicontrol('Parent',MC_panel,'units','normalized','position',[0.025,0.7,0.05,0.04],'style','text','string','0.0');
        North_string =uicontrol('Parent',MC_panel,'units','normalized','position',[0.45,0.85,0.05,0.04],'style','text','string','North');
        St_N=uicontrol('Parent',MC_panel,'units','normalized','position',[0.55,0.85,0.05,0.04],'style','text','String',num2str(get(S_N,'value')));

        S_S=uicontrol('Parent',MC_panel,'units','normalized','position',[0.05,0.55,0.9,0.05],'style','slider','Min',0.0,'Max',5.0,'value',P_S);
        Max_S=uicontrol('Parent',MC_panel,'units','normalized','position',[0.92,0.5,0.05,0.04],'style','text','string','5.0');
        Min_S=uicontrol('Parent',MC_panel,'units','normalized','position',[0.025,0.5,0.05,0.04],'style','text','string','0.0');
        South_string =uicontrol('Parent',MC_panel,'units','normalized','position',[0.45,0.65,0.05,0.04],'style','text','string','South');
        St_S=uicontrol('Parent',MC_panel,'units','normalized','position',[0.55,0.65,0.05,0.04],'style','text','String',num2str(get(S_S,'value')));

        S_E=uicontrol('Parent',MC_panel,'units','normalized','position',[0.05,0.35,0.9,0.05],'style','slider','Min',0.0,'Max',5.0,'value',P_E);
        Max_E=uicontrol('Parent',MC_panel,'units','normalized','position',[0.92,0.3,0.05,0.04],'style','text','string','5.0');
        Min_E=uicontrol('Parent',MC_panel,'units','normalized','position',[0.025,0.3,0.05,0.04],'style','text','string','0.0');
        East_string =uicontrol('Parent',MC_panel,'units','normalized','position',[0.45,0.45,0.05,0.04],'style','text','string','East');
        St_E=uicontrol('Parent',MC_panel,'units','normalized','position',[0.55,0.45,0.05,0.04],'style','text','String',num2str(get(S_E,'value')));

        S_W=uicontrol('Parent',MC_panel,'units','normalized','position',[0.05,0.15,0.9,0.05],'style','slider','Min',0.0,'Max',5.0,'value',P_E);
        Max_W=uicontrol('Parent',MC_panel,'units','normalized','position',[0.92,0.1,0.05,0.04],'style','text','string','5.0');
        Min_W=uicontrol('Parent',MC_panel,'units','normalized','position',[0.025,0.1,0.05,0.04],'style','text','string','0.0');
        West_string =uicontrol('Parent',MC_panel,'units','normalized','position',[0.45,0.25,0.05,0.04],'style','text','string','West');
        St_W=uicontrol('Parent',MC_panel,'units','normalized','position',[0.55,0.25,0.05,0.04],'style','text','String',num2str(get(S_W,'value')));
        
        IFT_panel=uipanel('Parent',h,'Title','IFT & Coalescence','background','white','FontSize',10,'Position',[0.03 0.42 0.4 0.08]);
        S_IFT=uicontrol('Parent',IFT_panel,'units','normalized','position',[0.05,0.45,0.9,0.25],'style','slider','Min',0.0,'Max',5.0,'value',P_N);
        Max_IFT=uicontrol('Parent',IFT_panel,'units','normalized','position',[0.92,0.1,0.05,0.25],'style','text','string','5.0');
        Min_IFT=uicontrol('Parent',IFT_panel,'units','normalized','position',[0.025,0.1,0.05,0.25],'style','text','string','0.0');
        St_IFT=uicontrol('Parent',IFT_panel,'units','normalized','position',[0.5,0.75,0.05,0.25],'style','text','String',num2str(get(S_IFT,'value')));
        
        N_string =uicontrol('Parent',h,'units','normalized','position',[0.95,0.90,0.015,0.015],'style','text','string','N');
        S_string =uicontrol('Parent',h,'units','normalized','position',[0.485,0.245,0.015,0.015],'style','text','string','S');
        E_string =uicontrol('Parent',h,'units','normalized','position',[0.95,0.245,0.015,0.015],'style','text','string','E');
        W_string =uicontrol('Parent',h,'units','normalized','position',[0.485,0.90,0.015,0.015],'style','text','string','W');

        h9=imshow(im);

%% start acquisition again

start(vid);
i_track=1;

%% Control loop starts
tic
while (1>0)
   
    if(get(h2,'value')==1) 
        if(get(Toggle_dropdown_box,'value')==1)
            set(src,'Exposure',str2double(get(Exposure,'String')));

            P_N = get(S_N,'value');set(St_N,'String',P_N);
            P_S = get(S_S,'value');set(St_S,'String',P_S);
            P_E = get(S_E,'value');set(St_E,'String',P_E);
            P_W = get(S_W,'value');set(St_W,'String',P_W);

            V1=P_S;
            V2=P_W;
            V3=P_N;
            V4=P_E;

            s.outputSingleScan([V1,V2,V3,V4]);

            Kp = str2double(get(Gain,'string'));

            flag=0;
            flag_i_track=0;

            image=(15*getsnapshot(vid));flushdata(vid);
            low_in = 799/65535;
            high_in = 6181/65535;
            image = imadjust(image,[low_in high_in]);
            imsub = image-imref;
            set(h9,'CData',imsub);
            drawnow

        end

        if(get(Toggle_dropdown_box,'value')==2)    

            Kp = str2double(get(Gain,'string')); 
            box = 100;

            if flag==0
                P_click = ginput(1);
                xp=P_click(1,1); yp=P_click(1,2); 
            end

            if flag_i_track==0
                i_track=1;
            end

            flag_i_track=1;
            flag=1;

            image=(15*getsnapshot(vid));flushdata(vid);
            low_in = 799/65535;
            high_in = 6181/65535;
            image = imadjust(image,[low_in high_in]);
            imsub = image-imref;
            set(h9,'CData',imsub);
            drawnow

            imblackwhite = im2bw(imsub,thresh);

            cropi = imblackwhite(yp-box:yp+box,xp-box:xp+box);
            neworigin = [yp-box xp-box];

            c=cell2mat((struct2cell(regionprops(cropi,'centroid','Area'))'));

            [m n] = size(c);

            if m==0
                flag=0;
                continue;

            elseif m>1
                centers_cropped =  c(1:m,2:3);
                centers_final_P1 = centers_cropped(:,1)+ neworigin(2);
                centers_final_P2= centers_cropped(:,2)+ neworigin(1);

                distance = [centers_final_P1-xp centers_final_P2-yp];
                distance = sum(distance'.^2)';
                ind = find(distance==min(distance));

                if length(ind)>1
                    ind = ind(1);
                end

                centers_final=[centers_final_P1 centers_final_P2];
                curr_pos_click = centers_final(ind,:);
                xp = curr_pos_click(1);
                yp = curr_pos_click(2);

            else
                centers_cropped =  c(1:m,2:3);
                centers_final_P1 = centers_cropped(:,1)+ neworigin(2);
                centers_final_P2= centers_cropped(:,2)+ neworigin(1);
                centers_final=[centers_final_P1 centers_final_P2];
                curr_pos_click = centers_final(1,:);

                xp = curr_pos_click(1);
                yp = curr_pos_click(2);

            end

     

            particle_pos_1 = (curr_pos_click) - center' + dim/2;
            particle_pos_1 = [particle_pos_1(1) dim-particle_pos_1(2)];
            
            hold off;
            h9=imshow(imsub);hold on;
            plot(curr_pos_click(1),curr_pos_click(2),'go');

            X=particle_pos_1; %% analytical soln frame of reference

            X_norm1=X/dim; %cartesian ;  %% analytical soln frame of reference

            X_dummy(:,:,i_track) = X_norm1;                           

            Pg = [V1 V2 V3 V4];           
            stag_pos=solxy/dim;

            Pg_new(i_track,:)=Pg;               
            initial_pos = X_norm1;  %cartesian ;  %% analytical soln frame of reference

            curr_pos(:,:,i_track)=initial_pos; %cylindrical

            P0 = sum(Pg)/4;
            Q = Pg-P0; 
            Q = Q(1:end-1);

            q(1,:)=[Q, -sum(Q)];

            curr_pos=initial_pos;  %% this line loses history of the particle

            % Revision note: the variable 'i' here can be used as a 
            % tracking variable, but currently is unused as it is reset 
            % in 'i=1'.

            i=1;
            i=i+1;

            distance= stag_pos-curr_pos;

            u=Kp*distance(1); 
            v=Kp*distance(2);

            parameter= [curr_pos(1,1,i-1) curr_pos(1,2,i-1) u v];

            psi_new = fsolve(@non_linear_eqns,[0 0],options,parameter);

            Q1=1;
            Q2=psi_new(1)-Q1;
            Q3=psi_new(2)-Q1-Q2;
            Q4=-(Q1+Q2+Q3);
            Q_calc= [Q1 Q2 Q3 Q4];

%             high_in = 6181/65535;
%             image = imadjust(image,[low_in high_in]);
%             imsub = image-imref;
%             set(h9,'CData',imsub);
%             drawnow

            P_new(i,:)= fminsearch(@eqns_for_flowratepressure,[ 1 1 1 1],options,Q_calc);


            while (max(P_new(i,:))>1|| min(P_new(i,:))<0)

                if max(P_new(i,:))>1
                    P_new(i,:)=P_new(i,:)/max(P_new(i,:))*1;
                end

                if min(P_new(i,:))<0
                    P_new(i,:)=P_new(i,:)-min(P_new(i,:));
                end  

            end

            V1=P_new(i,1);
            V2=P_new(i,2);
            V3=P_new(i,3);
            V4=P_new(i,4); 

            s.outputSingleScan([V1,V2,V3,V4]);

%             pause(0.8)

            P_C=(V1+V2+V3+V4)/4;
            Q1=(V1-P_C); 
            Q2=(V2-P_C);
            Q3 = (V3-P_C);
            Q4= (V4-P_C);

            Xs=(fminsearch(@uvvalues_al,[0.5,0.5],options,[Q1 Q1+Q2 Q1+Q2+Q3])-0.5)*dim; 

            Xs(1)=center(1)+Xs(1);
            Xs(2)=center(2)-Xs(2);

            hold on;
            plot(center(1),center(2),'m*',Xs(1),Xs(2),'g*')            

            X_new(:,:,i_track)=X_dummy(:,:,i_track)*dim-solxy;

            error(i_track) = sqrt(X_new(1,1,i_track).^2+X_new(1,2,i_track).^2);
toc
        end
    end
    
   if (get(h2,'value')==2);
       
       IFT_track = IFT_track + 1;
       
       set(src,'Exposure',str2double(get(Exposure,'String')));
       
       xs=center(1);
       ys=center(2);
       
       P_IFT = get(S_IFT,'value');
       set(St_IFT,'String',P_IFT);
       
       P_E=str2double(get(P_MIN,'String'));
       P_W=str2double(get(P_MIN,'String'));
       
%        if(get(Mode_input,'value')==1)
           image=(15*getsnapshot(vid));flushdata(vid);
%             low_in = 799/65535;
%             high_in = 6181/65535;
            image = imadjust(image,[low_in high_in]);
            imsub = image-imref;
            set(h9,'CData',imsub);
            drawnow
            
            
           [X xp yp]=extract_particle_position(xs,ys,low_in,high_in,thresh,imref,image,xp,yp,box);
           
           curr_pos_click = [xp yp];
           hold off;
           h9=imshow(imsub);hold on;
             particle_pos_1 = (curr_pos_click) - center' + dim/2;
            particle_pos_1 = [particle_pos_1(1) dim-particle_pos_1(2)];

            if(get(symbols_on,'value')==1)
                plot(curr_pos_click(1),curr_pos_click(2),'go');
            end
            
           
           
           x_r=X(1);y_r=X(2);
           K=str2double(get(GainIFT,'String'));
           V=corrected_voltages(x_r,y_r,K,P_E,P_W,P_IFT,P_IFT);
           
           V1=V(1);
           V2=V(2);
           V3=V(3);
           V4=V(4);
           
       
            
%         if (get(Mode_input,'value')==2)
%             X=extract_particle_position(xs,ys,i,thresh);
%             X_norm=([1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2)]*X)/dim+0.5;
%             K=str2double(get(Gain,'String'));
%             [V1,V2,V3,V4]=stag_point_control(xs,ys,X_norm,K,V1,V2,V3,V4,P_IFT,str2double(get(P_MIN,'String')))
%         end
        
        s.outputSingleScan([V1,V2,V3,V4]);
        
        
         P_C=(V1+V2+V3+V4)/4;
            Q1=(V1-P_C); 
            Q2=(V2-P_C);
            Q3 = (V3-P_C);
            Q4= (V4-P_C);
        
        Xs=(fminsearch(@uvvalues_al,[0.5,0.5],options,[Q1 Q1+Q2 Q1+Q2+Q3])-0.5)*dim; 

            Xs(1)=center(1)+Xs(1);
            Xs(2)=center(2)-Xs(2);

            hold on;
            if(get(symbols_on,'value')==1)
                plot(center(1),center(2),'m*',Xs(1),Xs(2),'g*')
            end

%             X_new(:,:,i_track)=X_dummy(:,:,i_track)*dim-solxy;

%             error(i_track) = sqrt(X_new(1,1,i_track).^2+X_new(1,2,i_track).^2);

if (get(ift_image_saving_on, 'value')==1)
    
    
    cropi = imsub(yp-box:yp+box,xp-box:xp+box);
P_string = num2str(V);

fid = fopen([sprintf('exp_1_image%d (',i_track) P_string ').bin'],'w');
fwrite(fid,cropi,'uint16');
fclose(fid);
end

Pressure_track(IFT_track,:) = V;
 
        
    end
    
    if (get(h2,'value')==6)
    
        disp('Program will now exit');
        s.outputSingleScan([2,2,0,0]);
        stop(vid);
        flushdata(vid);
        delete(vid);
        save data Voltage_data
        save manual_data data
        clear vid
        warning('on','all');
        break;
        
    end                
   
    try
        title(h1,['North = ' num2str(V1) ', South = ' num2str(V3) ', East = ' num2str(V4) ' , West = ' num2str(V2) ', Iteration # ', num2str(i_track) ', u ' num2str(u) ',v ' num2str(v)] )
    catch
        title(h1,['North = ' num2str(V1) ', South = ' num2str(V3) ', East = ' num2str(V4) ' , West = ' num2str(V2) ', Iteration # ', num2str(i_track)] )
    end

    i_track=i_track + 1;
    
end

