

dyn_info = [0:2*n+1; x_coord';y_coord';d_vec';q_vec';e_vec';l_vec'];

fileID = fopen('D:\VRP\Matlabstuff\codes\myAlgo-Dynamic\dyn_input.txt','w');
%fprintf(fileID,'%6s %12s\r\n','x','exp(x)');
fprintf(fileID,'%2i %6.2f %6.2f %2i %2i %6.2f %6.2f \r\n',dyn_info);
fclose(fileID);