function h = circle(x,y,r)
    % used for plotting a BLACK circle on top of a 2D diagram
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, 'k','LineWidth',2);
    hold off
end