function [particle_pos xp yp] = extract_particle_position(xs,ys,low_in,high_in,thresh,imref,image,xp,yp,box)
% global No_drop handle2 drop_pos
global vid
% drop_positions=get_drop_positions(i,thresh);
% 
% while(isempty(drop_positions))
%     if (get(No_drop,'value')==1), break, end
%     drop_positions=get_drop_positions(i,thresh);
% end
% 
% distance=(drop_positions(:,1)-xs).^2+(drop_positions(:,2)-ys).^2;
% index_min_distance=(distance==min(distance));
% desired_drop_position=drop_positions(index_min_distance,1:2);drop_pos=desired_drop_position;


            
            
             imsub = image-imref;

            imblackwhite = im2bw(imsub,thresh);

            cropi = imblackwhite(yp-box:yp+box,xp-box:xp+box);
            neworigin = [yp-box xp-box];

            c=cell2mat((struct2cell(regionprops(cropi,'centroid','Area'))'));

            [m n] = size(c);
            
      
     if m==0
                flag=0;
%                 continue;

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
            
desired_drop_position = [xp yp];


particle_pos=[1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)]*(desired_drop_position-[xs,ys])';
% try
%     delete(handle2);
% end
% hold on;handle2=plot(desired_drop_position(:,1),desired_drop_position(:,2),'y.'); hold off;
