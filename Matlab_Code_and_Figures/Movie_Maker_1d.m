function Movie_Maker_1d(movie_plot,X,no_of_plots,file_name)

title = strcat(file_name,'.avi');

writerObj = VideoWriter(title);
open(writerObj);

figure(1)

plot(X,squeeze(movie_plot(1,:)),'black')
axis([min(X) max(X) min(squeeze(movie_plot(no_of_plots-1,:))) max(squeeze(movie_plot(no_of_plots-1,:)))])

set(gca,'nextplot','replacechildren');

for mm = 1:no_of_plots-1
    
    plot(X,squeeze(movie_plot(mm,:)),'black')
    axis([min(X) max(X) min(squeeze(movie_plot(no_of_plots-1,:))) max(squeeze(movie_plot(no_of_plots-1,:)))])

    frame = getframe;
    writeVideo(writerObj,frame);

end

close(writerObj);
