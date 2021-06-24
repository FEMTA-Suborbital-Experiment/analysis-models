clear
clc
all_SPARTA=load('surfIDs2_notitles.out');
IDs=all_SPARTA(:,1);
xVertex1_SPARTA=all_SPARTA(:,2);
yVertex1_SPARTA=all_SPARTA(:,3);
zVertex1_SPARTA=all_SPARTA(:,4);
xVertex2_SPARTA=all_SPARTA(:,5);
yVertex2_SPARTA=all_SPARTA(:,6);
zVertex2_SPARTA=all_SPARTA(:,7);
xVertex3_SPARTA=all_SPARTA(:,8);
yVertex3_SPARTA=all_SPARTA(:,9);
zVertex3_SPARTA=all_SPARTA(:,10);

x_vertices=[xVertex1_SPARTA xVertex2_SPARTA xVertex3_SPARTA];
y_vertices=[yVertex1_SPARTA yVertex2_SPARTA yVertex3_SPARTA];
z_vertices=[zVertex1_SPARTA zVertex2_SPARTA zVertex3_SPARTA];
% xVertex_emitting=all_STL(:,2);
% yVertex_emitting=all_STL(:,3);
% zVertex_emitting=all_STL(:,4);
emitting_surf=zeros(length(xVertex1_SPARTA),1);
count=1;

plot3(x_vertices,y_vertices,z_vertices)

% for i=1:length(xVertex_SPARTA)
%     if xVertex_SPARTA(i)==xVertex_emitting(i)
%         if yVertex_SPARTA(i)==yVertex_emitting(i)
%             if zVertex_SPARTA(i)==zVertex_emitting(i)
%                 emitting_surf(count)=IDs(i);
%                 count=count+1;
%             end
%         end
%     end
% end


