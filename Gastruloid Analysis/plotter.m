function plotter(ri_CC, num_objects, nuc_info, border, drug_name, drug_str, col_name)
    L = labelmatrix(ri_CC);
    RGB2 = label2rgb(L,'spring','k','noshuffle'); 

    figure(2)
    subplot(1,3,1)

    imshow(RGB2)
    hold on
    for ii = 1:num_objects
        edge = nuc_info(ii).ClosestEdge;
        cent = nuc_info(ii).Centroid;
        hline = line([edge(1), cent(1)], [edge(2), cent(2)]); 
        hline.LineStyle = '-';
        hline.Color = 'c';
        plot(border(:,1),border(:,2),'w','LineWidth',3)
    end
    hold off

    ylim([0 1024])
    xlim([0 1024])
    set(gca,'XTick',[], 'YTick', [])


    subplot(1,3,2)
    hold on
    plot([sub_struct.(col_name).C2(:).Dist2Edge],[sub_struct.(col_name).C2(:).MeanIntensity],'o')
    plot([sub_struct.(col_name).C3(:).Dist2Edge],[sub_struct.(col_name).C3(:).MeanIntensity],'o')
    plot([sub_struct.(col_name).C4(:).Dist2Edge],[sub_struct.(col_name).C4(:).MeanIntensity],'o')
    xlabel('Distance to Edge')
    ylabel('Mean Intensity')
    legend(Cnames(chan_idx))
    hold off



    subplot(1,3,3)
    view(3)
    plot3([sub_struct.(col_name).C2(:).MeanIntensity],[sub_struct.(col_name).C3(:).MeanIntensity],[sub_struct.(col_name).C4(:).MeanIntensity],'o');
    xlabel('Mean Intensity C2')
    ylabel('Mean Intensity C3')
    zlabel("Mean Intensity C4")
    hold off

    sgtitle(sprintf(drug_name{i}))
    set(gcf, 'Position',  [100, 100, 1000, 300])
    savefig(gcf, drug_str{i});
    close all