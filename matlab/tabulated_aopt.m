%% Write a_opt on file

k = 30;
a_opt = gen_aopt(k);
fid = fopen("aopt.dat","w+");
for i=1:k
    fprintf(fid,"%1.16f\n",a_opt(i));
end
fclose(fid);