count=1;
X_d=X_dummy*dim;
X_du=[X_d(1,1,:)+center(1),-X_d(1,2,:)+center(2);X_d(2,1,:)+center(1),-X_d(2,2,:)+center(2)]; % converst back to cooridnates sytsem of image
 for i=1:1:407
% for i=726:1:726
    i
% i=90;
fid=fopen([sprintf('exp_1_image%d (',i_track) P_string ').bin']);
al=fread(fid,[1200 1600],'*uint16');
fclose(fid);
al1=(uint8(double(al)/65535*255));
al2=im2uint16(al1);

% al3 = al2+imref;
al3=al2;
low_in = 0/65535;
high_in = 38000/65535;
al4 = imadjust(al3,[low_in high_in]);
figure(1000);imshow(al4);
hold on
% figure(1000);plot(X_du(1,1,i),X_du(1,2,i),'ro','MarkerSize',10,X_du(2,1,i),X_du(2,2,i),'bo','MarkerSize',10,solxy(1,1)+center(1),solxy(1,2)+center(2),'rs',solxy(2,1)+center(1),solxy(2,2)+center(2),'bs');
figure(1000);plot(X_du(1,1,i),X_du(1,2,i),'yo','MarkerSize',12,'LineWidth',2)
figure(1000);plot(X_du(2,1,i),X_du(2,2,i),'go','MarkerSize',10,'LineWidth',2)
figure(1000);plot(solxy(1,1)+center(1),solxy(1,2)+center(2),'ys','MarkerSize',12,'LineWidth',2)
figure(1000); plot(solxy(2,1)+center(1)+10,solxy(2,2)+center(2)-5,'gs','MarkerSize',12,'LineWidth',2);
hold off
movie(count)=getframe(gcf);

img(:,:,:,count)=frame2im(movie(count));

count=count+1;

 end